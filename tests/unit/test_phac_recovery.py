import csv
import gzip
from pathlib import Path

import pytest

from phaatlas.pipeline import phac_recovery as pr


def _write_matrix(path, rows):
    fieldnames = ["genome", "n_phaC", "n_families_present"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def test_identify_qualifying_genomes_filters_correctly(tmp_path):
    path = tmp_path / "matrix.tsv"
    _write_matrix(path, [
        {"genome": "G_QUALIFIES", "n_phaC": "0", "n_families_present": "6"},
        {"genome": "G_HAS_PHAC", "n_phaC": "2", "n_families_present": "8"},
        {"genome": "G_TOO_FEW", "n_phaC": "0", "n_families_present": "5"},
    ])
    result = pr.identify_qualifying_genomes(path, min_other_genes=5)
    assert result == {"G_QUALIFIES"}


def _write_family_metadata(path, rows):
    fieldnames = ["target_id", "pha_family", "genome"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def test_collect_other_pha_hits_filters_to_qualifying_genomes(tmp_path):
    phaa_path = tmp_path / "phaA_unique_targets_with_metadata.tsv"
    _write_family_metadata(phaa_path, [
        {"target_id": "T1", "pha_family": "phaA", "genome": "G_QUALIFIES"},
        {"target_id": "T2", "pha_family": "phaA", "genome": "G_NOT_QUALIFYING"},
    ])
    phaz_path = tmp_path / "phaZ_unique_targets_with_metadata.tsv"
    _write_family_metadata(phaz_path, [
        {"target_id": "T3", "pha_family": "phaZ", "genome": "G_QUALIFIES"},
    ])

    hits = pr.collect_other_pha_hits({"G_QUALIFIES"}, {"phaA": phaa_path, "phaZ": phaz_path})
    assert len(hits) == 2
    assert {(h.family, h.target_id) for h in hits} == {("phaA", "T1"), ("phaZ", "T3")}


def test_collect_other_pha_hits_rejects_phac_as_anchor(tmp_path):
    with pytest.raises(AssertionError):
        pr.collect_other_pha_hits({"G"}, {"phaC": tmp_path / "nonexistent.tsv"})


def test_write_genome_download_manifest(tmp_path):
    catalog_path = tmp_path / "OMDBv2.0_data.tsv"
    with open(catalog_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["GENOME", "SAMPLE", "STUDY", "GENOME_FILE", "GENES_NT_FILE", "GENES_AA_FILE", "GENES_GFF_FILE", "ANTISMASH_FILE"])
        writer.writerow(["G_QUALIFIES", "S1", "ST1", "url_genome", "url_nt", "url_aa_1", "url_gff", "url_as"])
        writer.writerow(["G_OTHER", "S2", "ST1", "url_genome", "url_nt", "url_aa_2", "url_gff", "url_as"])

    out_path = tmp_path / "manifest.tsv"
    n = pr.write_genome_download_manifest({"G_QUALIFIES"}, catalog_path, out_path)
    assert n == 1
    rows = list(csv.DictReader(open(out_path, newline=""), delimiter="\t"))
    assert rows == [{"genome": "G_QUALIFIES", "genes_aa_url": "url_aa_1"}]


def test_parse_prodigal_faa_groups_by_scaffold_preserving_order():
    records = [
        ("G-scaffold_1_1", "AAA"),
        ("G-scaffold_1_2", "BBB"),
        ("G-scaffold_2_1", "CCC"),
        ("G-scaffold_1_3", "DDD"),
    ]
    result = pr.parse_prodigal_faa(records)
    assert result == {
        "G-scaffold_1": [("G-scaffold_1_1", "AAA"), ("G-scaffold_1_2", "BBB"), ("G-scaffold_1_3", "DDD")],
        "G-scaffold_2": [("G-scaffold_2_1", "CCC")],
    }


def test_find_anchor_gene_ids_matches_by_exact_sequence():
    genome_genes = {
        "G-scaffold_1": [("G-scaffold_1_1", "AAA"), ("G-scaffold_1_2", "PHAA_SEQ"), ("G-scaffold_1_3", "CCC")],
    }
    found = pr.find_anchor_gene_ids(genome_genes, anchor_sequences={"PHAA_SEQ", "SOME_OTHER_ANCHOR_NOT_PRESENT"})
    assert found == ["G-scaffold_1_2"]


def test_extract_neighborhood_windows_clips_to_scaffold_and_dedupes():
    # 25-gene scaffold, anchor at index 12 (0-based) -- window=3 should give indices 9..15
    genes = [(f"G-scaffold_1_{i+1}", f"SEQ{i}") for i in range(25)]
    genome_genes = {"G-scaffold_1": genes}
    windows = pr.extract_neighborhood_windows(genome_genes, anchor_gene_ids=["G-scaffold_1_13"], window=3)
    expected_ids = {f"G-scaffold_1_{i+1}" for i in range(9, 16)}
    assert set(windows) == expected_ids

    # short scaffold (5 genes), window=10 -- should just return the whole scaffold, no error
    short_genes = [(f"G-scaffold_2_{i+1}", f"SEQ{i}") for i in range(5)]
    genome_genes_short = {"G-scaffold_2": short_genes}
    windows_short = pr.extract_neighborhood_windows(genome_genes_short, anchor_gene_ids=["G-scaffold_2_3"], window=10)
    assert set(windows_short) == {f"G-scaffold_2_{i+1}" for i in range(5)}


def test_extract_neighborhood_windows_merges_overlapping_anchors_without_duplicates():
    genes = [(f"G-scaffold_1_{i+1}", f"SEQ{i}") for i in range(30)]
    genome_genes = {"G-scaffold_1": genes}
    windows = pr.extract_neighborhood_windows(genome_genes, anchor_gene_ids=["G-scaffold_1_5", "G-scaffold_1_8"], window=2)
    # anchor@5(0-idx4): 2..6(1-idx3..7) ; anchor@8(0-idx7): 5..9(1-idx6..10) -> union 3..11 (1-indexed ids)
    expected = {f"G-scaffold_1_{i}" for i in range(3, 11)}
    assert set(windows) == expected
    assert len(windows) == len(expected)  # no duplicate-inflated count


def test_extract_neighborhood_windows_skips_unknown_anchor_silently():
    genome_genes = {"G-scaffold_1": [("G-scaffold_1_1", "AAA")]}
    windows = pr.extract_neighborhood_windows(genome_genes, anchor_gene_ids=["G-scaffold_1_NOT_REAL"], window=5)
    assert windows == {}


def test_write_neighborhood_fasta_round_trip(tmp_path):
    neighborhoods = {"G-scaffold_1_1": "A" * 65, "G-scaffold_1_2": "B" * 10}
    out_path = tmp_path / "neighborhood.faa"
    n = pr.write_neighborhood_fasta(neighborhoods, out_path)
    assert n == 2
    records = pr.read_fasta(out_path)
    assert dict(records) == neighborhoods
    text = out_path.read_text()
    # 65-char sequence should wrap at 60 chars/line
    assert "A" * 60 + "\n" + "A" * 5 in text


def _make_genome_fasta_bytes(records):
    text = "".join(f">{h}\n{s}\n" for h, s in records)
    return gzip.compress(text.encode())  # matches what default_fetch_gzipped_fasta decompresses


def test_run_neighborhood_extraction_batch_happy_path(tmp_path):
    genome_a_records = [(f"A-scaffold_1_{i+1}", f"SEQ{i}") for i in range(20)]
    genome_a_records[9] = ("A-scaffold_1_10", "ANCHOR_SEQ_A")  # 0-indexed position 9
    genome_b_records = [(f"B-scaffold_1_{i+1}", f"OTHER{i}") for i in range(5)]

    fake_data = {
        "http://x/A.faa.gz": _make_genome_fasta_bytes(genome_a_records),
        "http://x/B.faa.gz": _make_genome_fasta_bytes(genome_b_records),  # no anchor present
    }

    def fake_fetch(url):
        return gzip.decompress(fake_data[url])

    results = pr.run_neighborhood_extraction_batch(
        genome_urls={"A": "http://x/A.faa.gz", "B": "http://x/B.faa.gz"},
        genome_anchor_sequences={"A": {"ANCHOR_SEQ_A"}, "B": {"SOME_SEQ_NOT_IN_B"}},
        out_path=tmp_path / "neighborhoods.faa",
        window=3,
        fetch=fake_fetch,
        retries=1,
        log=lambda msg: None,
    )
    by_genome = {r.genome: r for r in results}
    assert by_genome["A"].status == "ok"
    assert by_genome["A"].n_neighborhood_genes == 7  # window=3 around index 9 -> indices 6..12
    assert by_genome["B"].status == "no_anchors_found"

    written = pr.read_fasta(tmp_path / "neighborhoods.faa")
    assert len(written) == 7
    assert ("A-scaffold_1_10", "ANCHOR_SEQ_A") in written


def test_run_neighborhood_extraction_batch_retries_then_records_failure(tmp_path):
    calls = {"n": 0}

    def always_fails(url):
        calls["n"] += 1
        raise ConnectionError("simulated network failure")

    results = pr.run_neighborhood_extraction_batch(
        genome_urls={"A": "http://x/A.faa.gz"},
        genome_anchor_sequences={"A": {"ANCHOR"}},
        out_path=tmp_path / "neighborhoods.faa",
        fetch=always_fails,
        retries=3,
        retry_delay_s=0,
        log=lambda msg: None,
    )
    assert calls["n"] == 3
    assert results[0].status == "download_failed"


def test_run_neighborhood_extraction_batch_skips_genome_with_no_recorded_anchors(tmp_path):
    results = pr.run_neighborhood_extraction_batch(
        genome_urls={"A": "http://x/A.faa.gz"},
        genome_anchor_sequences={},  # A has no recorded anchor sequences at all
        out_path=tmp_path / "neighborhoods.faa",
        fetch=lambda url: (_ for _ in ()).throw(AssertionError("should not fetch when no anchors")),
        log=lambda msg: None,
    )
    assert results[0].status == "no_anchors_found"


def test_write_batch_report_round_trip(tmp_path):
    results = [
        pr.BatchResult(genome="A", status="ok", n_neighborhood_genes=7),
        pr.BatchResult(genome="B", status="download_failed"),
    ]
    out_path = tmp_path / "report.tsv"
    pr.write_batch_report(results, out_path)
    rows = list(csv.DictReader(open(out_path, newline=""), delimiter="\t"))
    assert rows[0] == {"genome": "A", "status": "ok", "n_neighborhood_genes": "7"}
    assert rows[1] == {"genome": "B", "status": "download_failed", "n_neighborhood_genes": "0"}


def test_load_anchor_sequences_by_genome(tmp_path):
    hits_path = tmp_path / "other_pha_hits.tsv"
    with open(hits_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["genome", "family", "target_id"])
        writer.writerow(["G1", "phaA", "T1"])
        writer.writerow(["G1", "phaZ", "T2"])
        writer.writerow(["G2", "phaA", "T3"])  # T3 has no known sequence -- should be skipped, not KeyError

    result = pr.load_anchor_sequences_by_genome(hits_path, target_id_sequences={"T1": "SEQ1", "T2": "SEQ2"})
    assert result == {"G1": {"SEQ1", "SEQ2"}}


@pytest.mark.skipif(
    __import__("shutil").which("hmmsearch") is None,
    reason="hmmsearch not on PATH (conda activate pha_phylo)",
)
def test_run_hmmsearch_and_parse_against_real_hmms(tmp_path):
    """Live check against the two real HMMs built for this task
    (phac_recovery/hmm/) -- confirms a known-true phaC scores as a hit
    and an unrelated random sequence does not, at the relaxed default
    threshold this pipeline actually uses."""
    hmm_path = Path(__file__).resolve().parent.parent.parent / "phac_recovery" / "hmm" / "phaC_custom.hmm"
    if not hmm_path.exists():
        pytest.skip("phaC_custom.hmm not present (run phac_recovery/hmm/build_custom_phac_hmm.sh)")

    # real reference phaC (UniProt A0A011T1J1, already confirmed live earlier
    # against both HMMs at bit-score 727/218) -- not a fabricated sequence,
    # so a failure here means the HMM/search wiring broke, not bad test data
    true_phac_seq = (
        "MPRTIDESVSSEAYRAIDQMREALSAHLTGGLSPASLALALIDWYIHLAAAPGKRLELIDKAVRKSARLSTYLAAAGVDPETPPCIEPLPGDYRFRAQAW"
        "SQQPFSSFAQSFLLAQQWWHNATHEVPGVTPHHEDVVSFMARQMLDVFSPSNNPFTNPEVLSKTIEAGGTNFVKGFQNWMEDAARAASGRPPVGTENFTPG"
        "EEVAVTPGKVVYRNHLIELIQYEPATDKVLAEPVLIVPAWIMKYYILDLSPHNSLVRYLVEKGHTVFCISWRNPTADDRNLTMDDYRKLGVMAALDAINAI"
        "VPDRKIHATGYCLGGTLLSIAAAAMAQHRDDRLASLTLFAAQTDFAEPGELALFIDHSQLHFLESMMWNRGYLSADQMAGAFQLLRSNDLIWSRIVRDYML"
        "GERTPMNDLMAWNADSTRMPYRMHAEYLKRLYLDNELATARFMVDGRPAALQNVRVPMFVVGTERDHVAPWQSVYKIHHLTDTDLVFVLTSGGHNAGIVSE"
        "PGHKGRRYRFAERRQHDPYLDPQEWVETASARDGSWWLEWSDWLLGRSTQERVAPPATGAADKRYAPLEDAPGTYVFQR"
    )
    unrelated_seq = "MSISNIGIYYFGAISLVLIGLYGVLVKQFGISAALSTLIAGVLVTFFTQPWQIWKKRRG"

    fasta_path = tmp_path / "test.faa"
    fasta_path.write_text(f">TRUE_PHAC\n{true_phac_seq}\n>UNRELATED\n{unrelated_seq}\n")
    domtbl_path = tmp_path / "out.domtbl"

    from phaatlas.pipeline import phac_recovery as pr2
    pr2.run_hmmsearch(hmm_path, fasta_path, domtbl_path, dom_e=10.0)
    assert domtbl_path.exists()
    hits = pr2.parse_hmmsearch_domtbl(domtbl_path, profile_name="phaC_custom")
    hit_gene_ids = {h.gene_id for h in hits}
    assert "TRUE_PHAC" in hit_gene_ids
    assert "UNRELATED" not in hit_gene_ids


REAL_GENOME_FAA = Path("/tmp/omdb_catalog/test.faa")


@pytest.mark.skipif(not REAL_GENOME_FAA.exists(), reason="real downloaded OMDB test genome not present at /tmp/omdb_catalog/test.faa")
def test_end_to_end_against_real_downloaded_genome():
    """Integration check against a real OMDB genome
    (ACIN21-1_SAMN05421563_MAG_00000023, one of the actual 13,522
    qualifying genomes) -- takes one real gene's own sequence as a stand-in
    anchor (equivalent to what a target_db-extracted anchor sequence would
    look like) and confirms the whole pipeline mechanically works against
    real prodigal-header data, not just synthetic fixtures."""
    records = pr.read_fasta(REAL_GENOME_FAA)
    assert len(records) == 3108  # confirmed live via `grep -c '^>' test.faa`

    genome_genes = pr.parse_prodigal_faa(records)
    # pick a real gene from partway into the file as a stand-in "known anchor" sequence
    anchor_gene_id, anchor_seq = records[500]
    found = pr.find_anchor_gene_ids(genome_genes, anchor_sequences={anchor_seq})
    assert anchor_gene_id in found

    windows = pr.extract_neighborhood_windows(genome_genes, anchor_gene_ids=found, window=10)
    assert anchor_gene_id in windows
    assert 1 <= len(windows) <= 21  # up to 2*10+1, less if near a scaffold edge
