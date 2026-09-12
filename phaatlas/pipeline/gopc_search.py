"""Per-family MMseqs2 protein-protein search of the NR95 PHA reference set
against GOPC. Does NOT classify hits as true PHA proteins -- this stage
only retains alignment statistics and summarizes unique GOPC hits per
family; conflict resolution across families and true/false classification
are deliberately out of scope here.

Pipeline: export_query_fastas -> build_gopc_target_db (once) ->
run_family_search (per family) -> combine_all_families (once, after every
family of interest has been searched).
"""

from __future__ import annotations

import csv
import shutil
import statistics
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

FORMAT_OUTPUT_COLUMNS = [
    "query", "target", "evalue", "bits", "pident", "alnlen",
    "qstart", "qend", "qlen", "tstart", "tend", "tlen", "qcov", "tcov",
]

UNIQUE_TARGET_COLUMNS = [
    "target_id", "pha_family", "best_query", "best_evalue", "best_bitscore",
    "best_pident", "best_qcov", "best_tcov", "n_reference_queries_matching",
]

EVALUE_THRESHOLDS = [1e-5, 1e-10, 1e-20]
IDENTITY_QUANTILE_LABELS = ["min", "p10", "p25", "median", "p75", "p90", "max"]


class MMseqsNotFoundError(RuntimeError):
    pass


def _require_mmseqs(mmseqs_bin: str) -> None:
    if shutil.which(mmseqs_bin) is None:
        raise MMseqsNotFoundError(
            f"'{mmseqs_bin}' not found on PATH. Run ./cluster/install_mmseqs2.sh, "
            f"or pass --mmseqs-bin with the full path."
        )


def export_query_fastas(db_path: Path, queries_dir: Path) -> dict[str, int]:
    """Writes one <family>.faa per pha_family under queries_dir, headers
    containing ONLY the stable protein_id (no other metadata) -- family/
    organism/etc. is joined back afterward from SQLite/CSV, per spec.
    Source: NR95 representatives (protein_id == cluster95_representative)
    with a non-null sequence. A family with no NR95 representatives yet
    (not clustered, or genuinely no sequenced members) is skipped, not
    written as an empty file.
    """
    import sqlite3

    queries_dir.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            """
            SELECT protein_id, pha_family, sequence
            FROM protein_master_export
            WHERE sequence IS NOT NULL AND sequence != ''
              AND cluster95_representative IS NOT NULL
              AND protein_id = cluster95_representative
            ORDER BY pha_family, protein_id
            """
        ).fetchall()
    finally:
        conn.close()

    by_family: dict[str, list[sqlite3.Row]] = {}
    for row in rows:
        by_family.setdefault(row["pha_family"], []).append(row)

    counts: dict[str, int] = {}
    for family_id, family_rows in by_family.items():
        out_path = queries_dir / f"{family_id}.faa"
        with open(out_path, "w") as f:
            for row in family_rows:
                f.write(f">{row['protein_id']}\n")
                seq = row["sequence"]
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
        counts[family_id] = len(family_rows)
    return counts


def build_gopc_target_db(
    gopc_faa_path: Path,
    target_db_path: Path,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
    shuffle: bool = False,
    gpu_compatible: bool = False,
) -> None:
    """`mmseqs createdb <gopc fasta(.gz)> <target_db_path>` -- run once.
    MMseqs2 reads gzipped FASTA directly, so GOPC.geneset.pep.fa.gz never
    needs to be decompressed to disk separately for this step.

    shuffle=False by default: `createdb`'s own --shuffle defaults to true
    (it randomizes sequence order up front for better load-balancing when
    the target is later split across threads/nodes during a search) --
    confirmed live, at GOPC's real scale (hundreds of millions of
    sequences) this pushed createdb's memory use past 64GB and got the job
    OOM-killed. Disabling it trades away some of that load-balancing
    benefit for actually being able to build the database; pass
    shuffle=True if you have enough memory headroom to afford both.

    gpu_compatible=False by default: pass True to build with
    `--createdb-mode 2` (mmseqs' own GPU-compatible storage layout),
    required before `run_family_search(..., gpu=True)` can search against
    this database. This is a storage-format choice, not a CUDA runtime
    requirement -- confirmed via `mmseqs createdb --help` that the flag
    exists on the plain CPU build too, so this build step itself does NOT
    need to run on a GPU node even when building a GPU-compatible
    database. A database built one way can't be searched the other way;
    rebuild (same command, flip this flag) if you need to switch.
    """
    _require_mmseqs(mmseqs_bin)
    target_db_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [mmseqs_bin, "createdb", str(gopc_faa_path), str(target_db_path), "--shuffle", "1" if shuffle else "0"]
    if gpu_compatible:
        cmd += ["--createdb-mode", "2"]
    if threads:
        cmd += ["--threads", str(threads)]
    subprocess.run(cmd, check=True)


@dataclass
class FamilySearchSummary:
    family_id: str
    n_reference_queries: int
    n_alignments: int = 0
    n_unique_gopc_targets: int = 0
    n_targets_evalue: dict[float, int] = field(default_factory=dict)
    n_targets_qcov_ge_0_5: int = 0
    n_targets_bidir_cov_ge_0_5: int = 0
    identity_quantiles: dict[str, float | None] = field(default_factory=dict)
    query_hit_cap_warning: bool = False
    queries_at_cap: list[str] = field(default_factory=list)
    max_seqs: int = 0


def run_family_search(
    family_id: str,
    query_fasta: Path,
    target_db_path: Path,
    results_dir: Path,
    tmp_dir: Path,
    sensitivity: float = 7.0,
    evalue: float = 1e-5,
    coverage: float = 0.0,
    max_seqs: int = 10000,
    cap_warning_fraction: float = 0.95,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
    split_memory_limit: str | None = None,
    gpu: bool = False,
) -> FamilySearchSummary:
    """Runs `mmseqs easy-search` for one family, then builds
    <family>_unique_targets.tsv and <family>_summary.tsv from
    <family>_hits.tsv in a single streaming pass (hits.tsv can be large --
    up to n_queries * max_seqs rows for a saturated family -- so this never
    loads the whole file into memory at once, only the small per-target/
    per-query aggregates).

    split_memory_limit: mmseqs' own --split-memory-limit (e.g. "50G").
    Strongly recommended for a CPU (gpu=False) run under SLURM/any
    cgroup-limited scheduler -- confirmed live that omitting it against a
    target this large ("Query database size: ... / Error: Prefilter died /
    Error: Search died") gets the prefilter step OOM-killed: without an
    explicit limit, MMseqs2 sizes its target-database split against the
    NODE's total physical memory, not the job's actual cgroup allocation,
    so on a shared node it can assume far more headroom than it's really
    been given. Set this comfortably below whatever --mem the job
    requested (see cluster/run_gopc_search.sbatch, which does this
    automatically from SLURM_MEM_PER_NODE).

    gpu: pass --gpu 1 to mmseqs (confirmed present on `search`/
    `easy-search` in mmseqs2 18-8cc5c's own --help; the GPU-accelerated
    prefilter is described in Kallenborn et al., "GPU-accelerated homology
    search with MMseqs2," bioRxiv 2024.11.13.623350). Requires: a
    GPU-enabled mmseqs binary (mmseqs-linux-gpu, not the plain CPU build --
    see cluster/install_mmseqs2.sh), a target database built with
    build_gopc_target_db(..., gpu_compatible=True), and an actual GPU
    available to the process (CUDA_VISIBLE_DEVICES). split_memory_limit is
    still passed through if given, but it's a CPU-prefilter-era setting --
    the GPU path's own memory behavior (which we have not been able to
    validate end-to-end, since building/testing this locally would need an
    actual CUDA GPU) may or may not use it the same way.
    """
    _require_mmseqs(mmseqs_bin)
    results_dir.mkdir(parents=True, exist_ok=True)
    family_tmp = tmp_dir / family_id
    family_tmp.mkdir(parents=True, exist_ok=True)
    hits_path = results_dir / f"{family_id}_hits.tsv"

    cmd = [
        mmseqs_bin, "easy-search",
        str(query_fasta), str(target_db_path), str(hits_path), str(family_tmp),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "-c", str(coverage),
        "--max-seqs", str(max_seqs),
        "--alignment-mode", "3",
        "--format-mode", "4",  # BLAST-TAB + column headers, so hits.tsv is self-describing
        "--format-output", ",".join(FORMAT_OUTPUT_COLUMNS),
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    if split_memory_limit:
        cmd += ["--split-memory-limit", split_memory_limit]
    if gpu:
        cmd += ["--gpu", "1"]
    subprocess.run(cmd, check=True)

    n_reference_queries = _count_fasta_records(query_fasta)
    summary = FamilySearchSummary(family_id=family_id, n_reference_queries=n_reference_queries, max_seqs=max_seqs)

    best_by_target: dict[str, dict] = {}
    queries_by_target: dict[str, set[str]] = {}
    query_hit_counts: dict[str, int] = {}
    identities: list[float] = []

    with open(hits_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            summary.n_alignments += 1
            q, t = row["query"], row["target"]
            ev, bits, pident = float(row["evalue"]), float(row["bits"]), float(row["pident"])
            qcov, tcov = float(row["qcov"]), float(row["tcov"])

            query_hit_counts[q] = query_hit_counts.get(q, 0) + 1
            queries_by_target.setdefault(t, set()).add(q)
            identities.append(pident)

            current = best_by_target.get(t)
            if current is None or (ev, -bits) < (current["evalue"], -current["bits"]):
                best_by_target[t] = {
                    "target_id": t, "pha_family": family_id, "best_query": q,
                    "evalue": ev, "bits": bits, "pident": pident, "qcov": qcov, "tcov": tcov,
                }

    # Cap warning: any query landing at/near max_seqs likely had its true
    # hit set truncated by the prefilter, not genuinely exhausted -- flag
    # for a rerun at a higher --max-seqs (see module docstring / spec).
    cap_threshold = cap_warning_fraction * max_seqs
    summary.queries_at_cap = sorted(q for q, n in query_hit_counts.items() if n >= cap_threshold)
    summary.query_hit_cap_warning = bool(summary.queries_at_cap)

    unique_targets_path = results_dir / f"{family_id}_unique_targets.tsv"
    with open(unique_targets_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(UNIQUE_TARGET_COLUMNS)
        for t in sorted(best_by_target):
            b = best_by_target[t]
            writer.writerow([
                b["target_id"], b["pha_family"], b["best_query"],
                b["evalue"], b["bits"], b["pident"], b["qcov"], b["tcov"],
                len(queries_by_target[t]),
            ])

    summary.n_unique_gopc_targets = len(best_by_target)
    for thr in EVALUE_THRESHOLDS:
        summary.n_targets_evalue[thr] = sum(1 for b in best_by_target.values() if b["evalue"] <= thr)
    summary.n_targets_qcov_ge_0_5 = sum(1 for b in best_by_target.values() if b["qcov"] >= 0.5)
    summary.n_targets_bidir_cov_ge_0_5 = sum(1 for b in best_by_target.values() if b["qcov"] >= 0.5 and b["tcov"] >= 0.5)
    summary.identity_quantiles = _identity_quantiles(identities)

    _write_family_summary(results_dir / f"{family_id}_summary.tsv", summary)
    return summary


def _count_fasta_records(fasta_path: Path) -> int:
    n = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def _identity_quantiles(identities: list[float]) -> dict[str, float | None]:
    if not identities:
        return {label: None for label in IDENTITY_QUANTILE_LABELS}
    s = sorted(identities)
    if len(s) == 1:
        return {label: s[0] for label in IDENTITY_QUANTILE_LABELS}
    # statistics.quantiles(n=100) gives the 1..99 percentile cut points;
    # index k corresponds to the k-th percentile (index 0 == p1).
    pcts = statistics.quantiles(s, n=100, method="inclusive")
    return {
        "min": s[0],
        "p10": pcts[9],
        "p25": pcts[24],
        "median": statistics.median(s),
        "p75": pcts[74],
        "p90": pcts[89],
        "max": s[-1],
    }


def _write_family_summary(path: Path, summary: FamilySearchSummary) -> None:
    row = {
        "family_id": summary.family_id,
        "n_reference_queries": summary.n_reference_queries,
        "n_alignments": summary.n_alignments,
        "n_unique_GOPC_targets": summary.n_unique_gopc_targets,
    }
    for thr in EVALUE_THRESHOLDS:
        row[f"n_targets_evalue_{thr:g}"] = summary.n_targets_evalue.get(thr, 0)
    row["n_targets_qcov_ge_0.5"] = summary.n_targets_qcov_ge_0_5
    row["n_targets_bidir_cov_ge_0.5"] = summary.n_targets_bidir_cov_ge_0_5
    for label in IDENTITY_QUANTILE_LABELS:
        key = "median_identity" if label == "median" else f"identity_{label}"
        row[key] = summary.identity_quantiles.get(label)
    row["max_seqs"] = summary.max_seqs
    row["query_hit_cap_warning"] = summary.query_hit_cap_warning
    row["queries_at_cap"] = "|".join(summary.queries_at_cap)

    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        for k, v in row.items():
            writer.writerow([k, v])


def combine_all_families(results_dir: Path) -> tuple[int, int]:
    """Concatenates every <family>_unique_targets.tsv into
    all_families_unique_targets.tsv, and every <family>_summary.tsv into
    all_families_summary.tsv (one row per family). Deliberately does NOT
    deduplicate a GOPC target hit by multiple families -- both rows are
    kept, per spec ("do not resolve proteins hit by multiple families yet").
    """
    unique_target_files = sorted(results_dir.glob("*_unique_targets.tsv"))
    combined_targets_path = results_dir / "all_families_unique_targets.tsv"
    n_target_rows = 0
    with open(combined_targets_path, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(UNIQUE_TARGET_COLUMNS)
        for path in unique_target_files:
            with open(path, newline="") as in_f:
                reader = csv.reader(in_f, delimiter="\t")
                next(reader, None)  # skip that file's own header
                for row in reader:
                    writer.writerow(row)
                    n_target_rows += 1

    summary_files = sorted(results_dir.glob("*_summary.tsv"))
    combined_summaries: list[dict[str, str]] = []
    all_metrics: list[str] = []
    for path in summary_files:
        with open(path, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)
            metrics = {}
            for metric, value in reader:
                metrics[metric] = value
                if metric not in all_metrics:
                    all_metrics.append(metric)
            combined_summaries.append(metrics)

    combined_summary_path = results_dir / "all_families_summary.tsv"
    with open(combined_summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(all_metrics)
        for metrics in combined_summaries:
            writer.writerow([metrics.get(m, "") for m in all_metrics])

    return n_target_rows, len(combined_summaries)
