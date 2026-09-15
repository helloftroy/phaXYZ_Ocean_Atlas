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

from phaatlas.db.session import connect_readonly

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

    Writes every family's file each call (not just one), and each write is
    write-to-temp-then-os.replace -- important because cluster/run_gopc_search.sbatch
    calls this at the start of EVERY per-family job, and multiple families'
    jobs commonly run concurrently on separate GPUs. Without the atomic
    swap, a concurrent mmseqs process reading e.g. queries/phaA.faa could
    observe a half-written file from an unrelated phaB job's own
    export-queries call; os.replace() on the same filesystem is atomic, so
    a reader always sees either the complete old file or the complete new
    one, never a partial write.
    """
    import os
    import sqlite3

    queries_dir.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
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
        tmp_path = queries_dir / f".{family_id}.faa.tmp.{os.getpid()}"
        with open(tmp_path, "w") as f:
            for row in family_rows:
                f.write(f">{row['protein_id']}\n")
                seq = row["sequence"]
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60] + "\n")
        os.replace(tmp_path, out_path)
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


def run_mmseqs_only(
    query_fasta: Path,
    target_db_path: Path,
    hits_path: Path,
    tmp_dir: Path,
    sensitivity: float = 7.0,
    evalue: float = 1e-5,
    coverage: float = 0.0,
    max_seqs: int = 10000,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
    split_memory_limit: str | None = None,
    gpu: bool = False,
    db_load_mode: int | None = None,
) -> None:
    """Runs `mmseqs easy-search` for one query FASTA against the GOPC
    target database, writing hits_path. No aggregation -- see
    aggregate_hits_file() for that, kept separate so a batched run
    (run_family_search_batch() below) can call this once per batch and
    aggregate only after every batch's hits file exists.

    db_load_mode: mmseqs' own --db-load-mode (0 auto / 1 fread / 2 mmap /
    3 mmap+touch), default None leaves mmseqs' own "auto" choice
    unchanged. Confirmed live: mmseqs has no dedicated memory cap for its
    final convertalis (output-formatting) step the way --split-memory-limit
    covers the prefilter, and a large family (phaA: 1822 queries, 18.2M
    alignments) was OOM-killed there even with a generous --mem. Explicit
    mmap (2) is the untested next thing to try if raising --mem alone
    isn't enough -- cgroups can reclaim clean mmap'd pages under memory
    pressure before resorting to an OOM-kill, unlike fread's fully-loaded
    buffers. Not defaulted on since it hasn't been verified to actually
    help this specific failure mode, only that it's a real, low-risk lever
    to reach for.

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

    IMPORTANT, confirmed live: mmseqs creates a fresh, RANDOMLY-named
    subdirectory under tmp_dir on every single invocation (not a
    deterministic one keyed by the query/target/params) -- resubmitting
    the exact same command after a job is killed does NOT resume from
    where the previous attempt left off, it restarts the entire scan of
    the target database from zero. There is no known way to force mmseqs
    to reuse a prior run's intermediate state through this CLI. If one
    run can't realistically finish inside your cluster's job time limit,
    use run_family_search_batch()/finalize_family_batches() below instead
    of resubmitting the same full-query command repeatedly.
    """
    _require_mmseqs(mmseqs_bin)
    hits_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        mmseqs_bin, "easy-search",
        str(query_fasta), str(target_db_path), str(hits_path), str(tmp_dir),
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
    if db_load_mode is not None:
        cmd += ["--db-load-mode", str(db_load_mode)]
    subprocess.run(cmd, check=True)


def aggregate_hits_file(
    hits_path: Path,
    family_id: str,
    n_reference_queries: int,
    results_dir: Path,
    max_seqs: int = 10000,
    cap_warning_fraction: float = 0.95,
) -> FamilySearchSummary:
    """Builds <family>_unique_targets.tsv and <family>_summary.tsv from an
    ALREADY-COMPLETE hits_path (either a single run's direct output, or
    several batches concatenated by finalize_family_batches()) in a single
    streaming pass -- hits.tsv can be large (up to n_queries * max_seqs
    rows for a saturated family), so this never loads the whole file into
    memory at once, only the small per-target/per-query aggregates.
    """
    results_dir.mkdir(parents=True, exist_ok=True)
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
    db_load_mode: int | None = None,
) -> FamilySearchSummary:
    """Single-shot search: run_mmseqs_only() against ALL of a family's
    queries in one mmseqs invocation, then aggregate_hits_file() on the
    result. Fine for a family whose full search comfortably fits inside
    one job's time limit; use run_family_search_batch()/
    finalize_family_batches() instead if it doesn't (see run_mmseqs_only's
    docstring for why resubmitting this doesn't help once a job is killed
    partway through)."""
    hits_path = results_dir / f"{family_id}_hits.tsv"
    run_mmseqs_only(
        query_fasta, target_db_path, hits_path, tmp_dir / family_id,
        sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
        mmseqs_bin=mmseqs_bin, threads=threads, split_memory_limit=split_memory_limit, gpu=gpu,
        db_load_mode=db_load_mode,
    )
    n_reference_queries = _count_fasta_records(query_fasta)
    return aggregate_hits_file(hits_path, family_id, n_reference_queries, results_dir, max_seqs, cap_warning_fraction)


def _read_fasta_records(fasta_path: Path) -> list[tuple[str, str]]:
    """Returns [(header_line_without_>, sequence), ...], sequence with
    internal newlines removed (single string, rewrapped on write)."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def _write_fasta_records(records: list[tuple[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def split_query_batch_fasta(query_fasta: Path, n_batches: int, batch_index: int, out_path: Path) -> int:
    """Writes out_path with every n_batches-th record starting at
    batch_index (simple round-robin striping -- keeps batches close to
    equal size regardless of any ordering in the source file, and is
    trivially deterministic/reproducible given the same n_batches). Returns
    the number of records written."""
    if not (0 <= batch_index < n_batches):
        raise ValueError(f"batch_index {batch_index} out of range for n_batches {n_batches}")
    records = _read_fasta_records(query_fasta)
    my_records = records[batch_index::n_batches]
    _write_fasta_records(my_records, out_path)
    return len(my_records)


def run_family_search_batch(
    family_id: str,
    query_fasta: Path,
    target_db_path: Path,
    results_dir: Path,
    tmp_dir: Path,
    n_batches: int,
    batch_index: int,
    sensitivity: float = 7.0,
    evalue: float = 1e-5,
    coverage: float = 0.0,
    max_seqs: int = 10000,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
    split_memory_limit: str | None = None,
    gpu: bool = False,
    db_load_mode: int | None = None,
) -> Path:
    """Runs mmseqs against ONLY this batch's slice of a family's queries
    (see split_query_batch_fasta), writing
    results_dir/<family_id>_batch<batch_index>of<n_batches>_hits.tsv.
    Deliberately does NOT aggregate -- run one call per batch_index
    (0..n_batches-1), normally as separate cluster jobs (e.g. a SLURM job
    array) each within its own time limit, THEN call
    finalize_family_batches() once every batch's hits file exists. Returns
    the hits file path written.

    This exists specifically because a single mmseqs invocation against
    the full GOPC target restarts from scratch if killed partway through
    (see run_mmseqs_only's docstring) -- splitting the QUERY side into
    batches doesn't reduce the cost of scanning the target once per batch
    (each batch still scans all of GOPC, so total GPU-time goes up, not
    down), but it turns "one run that might never finish inside the time
    limit" into several independent runs that each definitely can.
    """
    batch_dir = results_dir / "batches" / family_id
    batch_query_fasta = batch_dir / f"{family_id}_batch{batch_index}of{n_batches}.faa"
    n = split_query_batch_fasta(query_fasta, n_batches, batch_index, batch_query_fasta)

    hits_path = batch_dir / f"{family_id}_batch{batch_index}of{n_batches}_hits.tsv"
    batch_tmp = tmp_dir / family_id / f"batch{batch_index}of{n_batches}"
    run_mmseqs_only(
        batch_query_fasta, target_db_path, hits_path, batch_tmp,
        sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
        mmseqs_bin=mmseqs_bin, threads=threads, split_memory_limit=split_memory_limit, gpu=gpu,
        db_load_mode=db_load_mode,
    )
    print(f"{family_id} batch {batch_index}/{n_batches}: {n} queries -> {hits_path}")
    return hits_path


def finalize_family_batches(
    family_id: str,
    query_fasta: Path,
    results_dir: Path,
    n_batches: int,
    max_seqs: int = 10000,
    cap_warning_fraction: float = 0.95,
) -> FamilySearchSummary:
    """Once every batch 0..n_batches-1 has a completed hits file (see
    run_family_search_batch), concatenates them into the SAME
    <family>_hits.tsv / <family>_unique_targets.tsv / <family>_summary.tsv
    paths a single-shot run_family_search() would have produced -- output
    is indistinguishable in format from an unbatched run, just assembled
    from N independent mmseqs invocations instead of one. Raises
    FileNotFoundError naming whichever batch is missing/incomplete.
    """
    batch_dir = results_dir / "batches" / family_id
    batch_paths = [batch_dir / f"{family_id}_batch{i}of{n_batches}_hits.tsv" for i in range(n_batches)]
    missing = [p for p in batch_paths if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"{family_id}: {len(missing)}/{n_batches} batch hits file(s) missing -- run "
            f"run_family_search_batch for each first. Missing: {[str(p) for p in missing]}"
        )

    hits_path = results_dir / f"{family_id}_hits.tsv"
    hits_path.parent.mkdir(parents=True, exist_ok=True)
    with open(hits_path, "w", newline="") as out_f:
        header_written = False
        for batch_path in batch_paths:
            with open(batch_path, newline="") as in_f:
                header = in_f.readline()
                if not header_written:
                    out_f.write(header)
                    header_written = True
                for line in in_f:
                    out_f.write(line)

    n_reference_queries = _count_fasta_records(query_fasta)
    return aggregate_hits_file(hits_path, family_id, n_reference_queries, results_dir, max_seqs, cap_warning_fraction)


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
