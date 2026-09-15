import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Protein
from phaatlas.db.session import get_sessionmaker, init_db
from phaatlas.pipeline import gopc_search

MMSEQS_AVAILABLE = shutil.which("mmseqs") is not None


def test_identity_quantiles_empty():
    q = gopc_search._identity_quantiles([])
    assert all(v is None for v in q.values())


def test_identity_quantiles_single_value():
    q = gopc_search._identity_quantiles([42.0])
    assert all(v == 42.0 for v in q.values())


def test_identity_quantiles_monotonic_and_bounded():
    values = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
    q = gopc_search._identity_quantiles(values)
    assert q["min"] == 10.0
    assert q["max"] == 100.0
    assert q["median"] == 55.0
    ordered = [q["min"], q["p10"], q["p25"], q["median"], q["p75"], q["p90"], q["max"]]
    assert ordered == sorted(ordered)


def test_count_fasta_records(tmp_path):
    fasta = tmp_path / "q.faa"
    fasta.write_text(">a\nMKV\n>b\nMKL\n>c\nMKM\n")
    assert gopc_search._count_fasta_records(fasta) == 3


@pytest.fixture
def session(tmp_path):
    path = tmp_path / "test.sqlite"
    init_db(path)
    Session = get_sessionmaker(path)
    s = Session()
    s.add(FamilyDefinitionRow(family_id="phaC", display_name="PhaC", role="synthase", config_json="{}"))
    s.add(FamilyDefinitionRow(family_id="phaA", display_name="PhaA", role="thiolase", config_json="{}"))
    s.commit()
    yield s, path
    s.close()


def test_export_query_fastas_headers_are_protein_id_only(session):
    s, db_path = session
    s.add(Protein(protein_id="UNIPROT:X1", uniprot_accession="X1", sequence="MKVLLLL", sequence_sha256="h1"))
    s.flush()
    s.add(
        FamilyAssignment(
            protein_id="UNIPROT:X1", family_id="phaC", retrieval_reasons="uniprot_reviewed_gene_name",
            evidence_tier="GOLD_EXPERIMENTAL", needs_manual_review=False,
            cluster95_id="phaC:UNIPROT:X1", cluster95_representative="UNIPROT:X1", cluster95_size=1,
        )
    )
    s.commit()

    queries_dir = db_path.parent / "queries"
    counts = gopc_search.export_query_fastas(db_path, queries_dir)
    assert counts == {"phaC": 1}
    text = (queries_dir / "phaC.faa").read_text()
    assert text.splitlines()[0] == ">UNIPROT:X1"  # ONLY the protein_id, no family=/organism=/etc


def test_export_query_fastas_skips_non_representative_members(session):
    s, db_path = session
    s.add(Protein(protein_id="UNIPROT:R1", uniprot_accession="R1", sequence="MKVLLLL", sequence_sha256="h1"))
    s.add(Protein(protein_id="UNIPROT:M1", uniprot_accession="M1", sequence="MKVLLLL", sequence_sha256="h1"))
    s.flush()
    # R1 is the cluster representative; M1 shares its cluster but is NOT
    # the representative -- only R1 should appear in the query FASTA.
    for pid in ("UNIPROT:R1", "UNIPROT:M1"):
        s.add(
            FamilyAssignment(
                protein_id=pid, family_id="phaC", retrieval_reasons="uniprot_reviewed_gene_name",
                evidence_tier="GOLD_EXPERIMENTAL", needs_manual_review=False,
                cluster95_id="phaC:UNIPROT:R1", cluster95_representative="UNIPROT:R1", cluster95_size=2,
            )
        )
    s.commit()

    queries_dir = db_path.parent / "queries"
    counts = gopc_search.export_query_fastas(db_path, queries_dir)
    assert counts == {"phaC": 1}


def test_combine_all_families(tmp_path):
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    (results_dir / "phaA_unique_targets.tsv").write_text(
        "\t".join(gopc_search.UNIQUE_TARGET_COLUMNS) + "\n"
        "T1\tphaA\tQ1\t1e-10\t100\t80\t0.9\t0.9\t3\n"
    )
    (results_dir / "phaC_unique_targets.tsv").write_text(
        "\t".join(gopc_search.UNIQUE_TARGET_COLUMNS) + "\n"
        "T1\tphaC\tQ2\t1e-20\t200\t95\t1.0\t1.0\t5\n"  # SAME target T1 hit by a different family -- must NOT be deduped
        "T2\tphaC\tQ3\t1e-8\t90\t70\t0.6\t0.6\t1\n"
    )
    (results_dir / "phaA_summary.tsv").write_text("metric\tvalue\nfamily_id\tphaA\nn_alignments\t10\n")
    (results_dir / "phaC_summary.tsv").write_text("metric\tvalue\nfamily_id\tphaC\nn_alignments\t20\n")

    n_targets, n_families = gopc_search.combine_all_families(results_dir)
    assert n_targets == 3  # T1 appears twice (once per family) -- not collapsed
    assert n_families == 2

    combined_text = (results_dir / "all_families_unique_targets.tsv").read_text()
    assert combined_text.count("T1\t") == 2

    summary_text = (results_dir / "all_families_summary.tsv").read_text()
    assert "phaA" in summary_text and "phaC" in summary_text


@pytest.mark.skipif(not MMSEQS_AVAILABLE, reason="mmseqs binary not on PATH")
def test_run_family_search_end_to_end(tmp_path):
    query_fasta = tmp_path / "phaX.faa"
    seq = "MKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAFAAAAQPMDSPAALQAYTAKLGLPPAQAWTAQNFLES" * 2
    query_fasta.write_text(f">QUERY1\n{seq}\n")

    target_fasta = tmp_path / "target.faa"
    # a close match (few mismatches) and a clearly unrelated decoy
    close = seq[:-4] + "AAAA"
    decoy = "GGGGPPPPLLLLSSSSTTTTNNNNQQQQEEEEDDDDKKKKRRRRHHHHYYYYWWWWFFFF" * 5
    target_fasta.write_text(f">TARGET_CLOSE\n{close}\n>TARGET_DECOY\n{decoy}\n")

    target_db = tmp_path / "db" / "target_db"
    gopc_search.build_gopc_target_db(target_fasta, target_db)

    results_dir = tmp_path / "results"
    tmp_dir = tmp_path / "tmp"
    summary = gopc_search.run_family_search(
        "phaX", query_fasta, target_db, results_dir, tmp_dir,
        sensitivity=7.0, evalue=1e-5, coverage=0.0, max_seqs=10000,
    )

    assert summary.n_reference_queries == 1
    assert summary.n_unique_gopc_targets == 1  # only the close match passes e-1e-5, decoy shouldn't
    assert summary.query_hit_cap_warning is False

    unique_targets_text = (results_dir / "phaX_unique_targets.tsv").read_text()
    assert "TARGET_CLOSE" in unique_targets_text
    assert "TARGET_DECOY" not in unique_targets_text

    summary_text = (results_dir / "phaX_summary.tsv").read_text()
    assert "n_unique_GOPC_targets\t1" in summary_text


@pytest.mark.skipif(not MMSEQS_AVAILABLE, reason="mmseqs binary not on PATH")
def test_run_mmseqs_only_resumes_convertalis_without_rerunning_search(tmp_path, monkeypatch):
    """The whole point of the split search/convertalis refactor: if
    convertalis fails (e.g. OOM, as happened to phaA in production) after
    search already completed, a second call must reuse the existing
    alignment results rather than rescanning the target database again."""
    query_fasta = tmp_path / "phaX.faa"
    seq = "MKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAFAAAAQPMDSPAALQAYTAKLGLPPAQAWTAQNFLES" * 2
    query_fasta.write_text(f">QUERY1\n{seq}\n")

    target_fasta = tmp_path / "target.faa"
    close = seq[:-4] + "AAAA"
    decoy = "GGGGPPPPLLLLSSSSTTTTNNNNQQQQEEEEDDDDKKKKRRRRHHHHYYYYWWWWFFFF" * 5
    target_fasta.write_text(f">TARGET_CLOSE\n{close}\n>TARGET_DECOY\n{decoy}\n")

    target_db = tmp_path / "db" / "target_db"
    gopc_search.build_gopc_target_db(target_fasta, target_db)

    hits_path = tmp_path / "hits.tsv"
    tmp_dir = tmp_path / "tmp"

    gopc_search.run_mmseqs_only(query_fasta, target_db, hits_path, tmp_dir)
    assert (tmp_dir / "result.dbtype").exists()
    first_hits = hits_path.read_text()

    # Simulate convertalis having failed/been deleted after a successful
    # search: remove only the final TSV, not the alignment result db.
    hits_path.unlink()

    calls = []
    real_run = subprocess.run

    def spy_run(cmd, *args, **kwargs):
        calls.append(cmd)
        return real_run(cmd, *args, **kwargs)

    monkeypatch.setattr(gopc_search.subprocess, "run", spy_run)
    gopc_search.run_mmseqs_only(query_fasta, target_db, hits_path, tmp_dir)

    assert not any(cmd[1] == "search" for cmd in calls), "resumed run must not re-invoke mmseqs search"
    assert any(cmd[1] == "convertalis" for cmd in calls)
    assert hits_path.read_text() == first_hits


def test_split_query_batch_fasta_round_robin(tmp_path):
    query_fasta = tmp_path / "q.faa"
    query_fasta.write_text("".join(f">Q{i}\nMK{i}VL\n" for i in range(7)))

    out0 = tmp_path / "batch0.faa"
    out1 = tmp_path / "batch1.faa"
    out2 = tmp_path / "batch2.faa"
    n0 = gopc_search.split_query_batch_fasta(query_fasta, 3, 0, out0)
    n1 = gopc_search.split_query_batch_fasta(query_fasta, 3, 1, out1)
    n2 = gopc_search.split_query_batch_fasta(query_fasta, 3, 2, out2)

    assert (n0, n1, n2) == (3, 2, 2)  # 7 records over 3 batches: sizes as even as round-robin allows
    assert n0 + n1 + n2 == 7

    headers = lambda p: [l[1:].strip() for l in p.read_text().splitlines() if l.startswith(">")]
    all_headers = headers(out0) + headers(out1) + headers(out2)
    assert sorted(all_headers) == sorted(f"Q{i}" for i in range(7))  # every record assigned exactly once
    assert headers(out0) == ["Q0", "Q3", "Q6"]  # stride-3 starting at index 0


def test_split_query_batch_fasta_rejects_bad_index(tmp_path):
    query_fasta = tmp_path / "q.faa"
    query_fasta.write_text(">Q0\nMKV\n")
    with pytest.raises(ValueError):
        gopc_search.split_query_batch_fasta(query_fasta, 3, 3, tmp_path / "out.faa")


@pytest.mark.skipif(not MMSEQS_AVAILABLE, reason="mmseqs binary not on PATH")
def test_batched_search_matches_unbatched(tmp_path):
    """The whole point of batching: splitting queries across N separate
    mmseqs invocations and finalizing must produce the SAME
    unique_targets/summary as running all queries in one shot -- batching
    is purely how the compute is scheduled, not a change in what's found.
    """
    base_seq = "MKVLNRQAVASLKELQASAAAINSNPFAAAKPAEIQGLARFVQAAKADPAGAFAAAAQPMDSPAALQAYTAKLGLPPAQAWTAQNFLES"

    def mutate(seq: str, positions: list[int]) -> str:
        chars = list(seq)
        for p in positions:
            chars[p] = "A" if chars[p] != "A" else "G"
        return "".join(chars)

    queries = [(f"QUERY{i}", mutate(base_seq, [i, i + 10])) for i in range(5)]
    query_fasta = tmp_path / "queries.faa"
    query_fasta.write_text("".join(f">{h}\n{s}\n" for h, s in queries))

    targets = [(f"CLOSE{i}", mutate(base_seq, [i, i + 20])) for i in range(5)]
    targets.append(("DECOY", "GGGGPPPPLLLLSSSSTTTTNNNNQQQQEEEEDDDDKKKKRRRRHHHHYYYYWWWWFFFF" * 5))
    target_fasta = tmp_path / "targets.faa"
    target_fasta.write_text("".join(f">{h}\n{s}\n" for h, s in targets))

    target_db = tmp_path / "db" / "target_db"
    gopc_search.build_gopc_target_db(target_fasta, target_db)

    # Unbatched reference run.
    ref_results = tmp_path / "ref_results"
    ref_summary = gopc_search.run_family_search(
        "phaX", query_fasta, target_db, ref_results, tmp_path / "ref_tmp", max_seqs=10000,
    )

    # Batched run: 3 batches, each its own mmseqs invocation, then finalize.
    batch_results = tmp_path / "batch_results"
    for i in range(3):
        gopc_search.run_family_search_batch(
            "phaX", query_fasta, target_db, batch_results, tmp_path / "batch_tmp",
            n_batches=3, batch_index=i, max_seqs=10000,
        )
    batched_summary = gopc_search.finalize_family_batches("phaX", query_fasta, batch_results, n_batches=3)

    assert batched_summary.n_reference_queries == ref_summary.n_reference_queries
    assert batched_summary.n_alignments == ref_summary.n_alignments
    assert batched_summary.n_unique_gopc_targets == ref_summary.n_unique_gopc_targets

    def read_unique_targets(results_dir: Path) -> set[tuple]:
        with open(results_dir / "phaX_unique_targets.tsv", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            return {(r["target_id"], r["best_query"], r["n_reference_queries_matching"]) for r in reader}

    assert read_unique_targets(batch_results) == read_unique_targets(ref_results)
