from __future__ import annotations

import json
from pathlib import Path

import typer
from dotenv import load_dotenv
from rich.console import Console
from rich.table import Table
from sqlalchemy import func, select

from phaatlas.config_loader import DEFAULT_FAMILY_CONFIG_PATH, REPO_ROOT, load_family_definitions
from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Phenotype, Protein, RetrievalRun, SourceEvidence
from phaatlas.db.session import get_db_path, init_db, session_scope
from phaatlas.pipeline import export as export_pipeline
from phaatlas.pipeline import gopc_search as gopc_search_pipeline
from phaatlas.pipeline.cluster95 import MMseqsNotFoundError, cluster_family
from phaatlas.pipeline.ingest_brenda import ingest_brenda_for_family
from phaatlas.pipeline.ingest_uniprot import ingest_phaR_disambiguation, ingest_uniprot_for_family
from phaatlas.sources import brenda as brenda_source

load_dotenv()

app = typer.Typer(help="PHA reference database pipeline (BRENDA + UniProt, family-config driven).")
console = Console()


def _sync_family_definitions(session) -> None:
    for fam in load_family_definitions():
        row = session.get(FamilyDefinitionRow, fam.family_id)
        if row is None:
            row = FamilyDefinitionRow(family_id=fam.family_id)
            session.add(row)
        row.display_name = fam.display_name
        row.role = fam.role
        row.ec_brenda = fam.ec_brenda
        row.require_pha_context = fam.require_pha_context
        row.config_json = fam.raw_json()
    session.commit()


def _resolve_families(family_arg: str):
    all_families = load_family_definitions()
    if family_arg == "all":
        return all_families
    by_id = {f.family_id: f for f in all_families}
    requested = [f.strip() for f in family_arg.split(",") if f.strip()]
    unknown = [f for f in requested if f not in by_id]
    if unknown:
        console.print(f"[red]Unknown family_id(s) {unknown}. Known: {', '.join(sorted(by_id))}[/red]")
        raise typer.Exit(code=2)
    return [by_id[f] for f in requested]


@app.command("init-db")
def init_db_cmd():
    """Create/refresh all tables and the protein_master_export view, and
    sync family_definition from config/family_definitions.yaml."""
    init_db()
    with session_scope() as session:
        _sync_family_definitions(session)
    console.print(f"[green]Database ready at {get_db_path()}[/green]")


@app.command("ingest-uniprot")
def ingest_uniprot_cmd(
    family: str = typer.Option("phaC", help="family_id from config/family_definitions.yaml, or 'all'"),
    max_results: int = typer.Option(1000, help="max results per individual UniProt query"),
):
    """Route B: family-specific UniProt reviewed/unreviewed searches."""
    families = _resolve_families(family)
    phaR_ids = {"phaR_synthase", "phaR_regulator"}
    families_by_id = {f.family_id: f for f in families}

    # phaR_synthase and phaR_regulator share an exact gene name (phaR) and
    # MUST always be disambiguated jointly (see ingest_phaR_disambiguation) --
    # requesting only one of them still pulls in its partner here so a
    # `gene:phaR` hit is never blindly assigned to just the one requested,
    # which would silently skip the synthase/regulator disambiguation check.
    if families_by_id.keys() & phaR_ids and not phaR_ids.issubset(families_by_id.keys()):
        missing = phaR_ids - families_by_id.keys()
        console.print(f"[dim]Pulling in {missing} for joint phaR disambiguation.[/dim]")
        all_by_id = {f.family_id: f for f in load_family_definitions()}
        for m in missing:
            families.append(all_by_id[m])
            families_by_id[m] = all_by_id[m]

    with session_scope() as session:
        _sync_family_definitions(session)

        if phaR_ids.issubset(families_by_id.keys()):
            console.print("Running joint phaR disambiguation (gene:phaR)...")
            summary = ingest_phaR_disambiguation(
                session, families_by_id["phaR_synthase"], families_by_id["phaR_regulator"], max_results_per_query=max_results
            )
            console.print(json.dumps(summary))

        for fam in families:
            skip_gene = fam.family_id in phaR_ids and phaR_ids.issubset(families_by_id.keys())
            console.print(f"[bold]{fam.family_id}[/bold]: querying UniProt...")
            summary = ingest_uniprot_for_family(session, fam, max_results_per_query=max_results, skip_gene_queries=skip_gene)
            console.print(json.dumps(summary))


@app.command("download-brenda")
def download_brenda_cmd(force: bool = typer.Option(False, help="re-download even if a cached copy exists")):
    """Downloads the full BRENDA JSON dump (accepting the CC BY 4.0 license
    the same way the website's own download form does)."""
    path = brenda_source.download_json(force=force)
    release = brenda_source.get_release(path)
    console.print(f"[green]BRENDA JSON ready at {path} (release={release})[/green]")


@app.command("ingest-brenda")
def ingest_brenda_cmd(
    family: str = typer.Option("phaC", help="family_id from config/family_definitions.yaml, or 'all'"),
    brenda_json: Path = typer.Option(None, help="path to an already-extracted brenda_<release>.json; downloads if omitted"),
):
    """EC-anchored BRENDA ingestion for families with ec_brenda set (skips
    families with ec_brenda: null, per spec)."""
    json_path = brenda_json or brenda_source.download_json()
    release = brenda_source.get_release(json_path)
    console.print(f"Using BRENDA {json_path} (release={release})")

    families = [f for f in _resolve_families(family) if f.ec_brenda]
    if not families:
        console.print("[yellow]No families with ec_brenda set in this selection -- nothing to do.[/yellow]")
        return

    with session_scope() as session:
        _sync_family_definitions(session)
        for fam in families:
            console.print(f"[bold]{fam.family_id}[/bold] (EC {fam.ec_brenda}): streaming BRENDA record...")
            summary = ingest_brenda_for_family(session, fam, json_path, release)
            console.print(json.dumps(summary))


@app.command("export")
def export_cmd():
    """Regenerates exports/pha_reference_master.csv and pha_reference.faa
    from the SQLite database."""
    db_path = get_db_path()
    csv_path = db_path.parent / "exports" / "pha_reference_master.csv"
    faa_path = db_path.parent / "exports" / "pha_reference.faa"
    n_rows = export_pipeline.export_master_csv(db_path, csv_path)
    n_seqs = export_pipeline.export_fasta(db_path, faa_path)
    console.print(f"[green]{csv_path}: {n_rows} rows[/green]")
    console.print(f"[green]{faa_path}: {n_seqs} sequences[/green]")


@app.command("cluster95")
def cluster95_cmd(
    family: str = typer.Option("all", help="family_id, comma-separated list, or 'all'"),
    min_seq_id: float = typer.Option(0.95, help="MMseqs2 --min-seq-id"),
    min_cov: float = typer.Option(0.90, help="MMseqs2 -c (coverage, --cov-mode 0 i.e. bidirectional)"),
    mmseqs_bin: str = typer.Option("mmseqs", help="mmseqs binary name or full path"),
    threads: int = typer.Option(None, help="MMseqs2 --threads (default: mmseqs' own default)"),
):
    """Exact-dedups sequences then clusters WITHIN each family at
    min_seq_id/min_cov with MMseqs2, writing cluster95_id/
    cluster95_representative/cluster95_size back onto family_assignment.
    Never deletes rows -- only ADD COLUMN migrations + in-place UPDATEs.
    Run `export` afterward to regenerate the master CSV with these columns,
    or `export-cluster95` for the all/nr95/cluster-tsv bundle.
    """
    families = _resolve_families(family)
    with session_scope() as session:
        for fam in families:
            try:
                summary = cluster_family(
                    session, fam.family_id, min_seq_id=min_seq_id, min_cov=min_cov, mmseqs_bin=mmseqs_bin, threads=threads
                )
            except MMseqsNotFoundError as exc:
                console.print(f"[red]{exc}[/red]")
                raise typer.Exit(code=2)
            console.print(json.dumps(summary))


@app.command("export-cluster95")
def export_cluster95_cmd():
    """Writes exports/pha_reference_all.faa, pha_reference_nr95.faa, and
    pha_reference_cluster95.tsv from the SQLite database. Run `cluster95`
    first -- pha_reference_nr95.faa is empty for any family that hasn't
    been clustered yet."""
    db_path = get_db_path()
    exports_dir = db_path.parent / "exports"
    n_all = export_pipeline.export_fasta_all(db_path, exports_dir / "pha_reference_all.faa")
    n_nr95 = export_pipeline.export_fasta_nr95(db_path, exports_dir / "pha_reference_nr95.faa")
    n_tsv = export_pipeline.export_cluster95_tsv(db_path, exports_dir / "pha_reference_cluster95.tsv")
    console.print(f"[green]{exports_dir / 'pha_reference_all.faa'}: {n_all} sequences[/green]")
    console.print(f"[green]{exports_dir / 'pha_reference_nr95.faa'}: {n_nr95} sequences (cluster representatives)[/green]")
    console.print(f"[green]{exports_dir / 'pha_reference_cluster95.tsv'}: {n_tsv} rows[/green]")


GOPC_SEARCH_DIR = REPO_ROOT / "PHA_bioprospecting" / "gopc_search"


@app.command("export-queries")
def export_queries_cmd(
    queries_dir: Path = typer.Option(GOPC_SEARCH_DIR / "queries", help="output directory, one <family>.faa per family"),
):
    """Writes one query FASTA per pha_family from the NR95 cluster
    representatives (cluster95_representative rows) -- headers contain
    ONLY the protein_id, no other metadata. Run `cluster95`/`export-cluster95`
    first; a family with no NR95 representatives yet is skipped."""
    db_path = get_db_path()
    counts = gopc_search_pipeline.export_query_fastas(db_path, queries_dir)
    for family_id, n in sorted(counts.items()):
        console.print(f"[green]{queries_dir / f'{family_id}.faa'}: {n} queries[/green]")
    console.print(f"[bold]{sum(counts.values())} total NR95 query sequences across {len(counts)} families[/bold]")


@app.command("gopc-build-db")
def gopc_build_db_cmd(
    gopc_faa: Path = typer.Argument(..., help="GOPC FASTA, .gz is fine -- mmseqs reads gzip directly, no need to decompress"),
    target_db: Path = typer.Option(GOPC_SEARCH_DIR / "target_db" / "gopc_db", help="mmseqs target database path (prefix)"),
    mmseqs_bin: str = typer.Option("mmseqs", help="mmseqs binary name or full path"),
    threads: int = typer.Option(None, help="mmseqs --threads (default: mmseqs' own default)"),
    shuffle: bool = typer.Option(
        False,
        help="mmseqs createdb --shuffle -- defaults OFF: at GOPC's real scale this pushed memory past 64GB and got "
        "the job OOM-killed (confirmed live). Only turn on if you have generous memory headroom to spare for "
        "better target-split load-balancing during later searches.",
    ),
    gpu_compatible: bool = typer.Option(
        False,
        help="mmseqs createdb --createdb-mode 2 -- required before `gopc-search --gpu` can search this database. "
        "A storage-format choice, not a CUDA requirement to build (this step doesn't need an actual GPU), but a "
        "database built one way can't be searched the other way -- rebuild (flip this flag) to switch.",
    ),
):
    """`mmseqs createdb` on GOPC -- run once. Slow/IO-heavy given GOPC's
    real size (~184GB compressed); run via cluster/run_gopc_build_db.sbatch,
    not interactively."""
    try:
        gopc_search_pipeline.build_gopc_target_db(
            gopc_faa, target_db, mmseqs_bin=mmseqs_bin, threads=threads, shuffle=shuffle, gpu_compatible=gpu_compatible
        )
    except gopc_search_pipeline.MMseqsNotFoundError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=2)
    console.print(f"[green]GOPC target database ready at {target_db}[/green]")


@app.command("gopc-search")
def gopc_search_cmd(
    family: str = typer.Option("all", help="family_id, comma-separated list, or 'all' (searches every queries/<family>.faa present)"),
    queries_dir: Path = typer.Option(GOPC_SEARCH_DIR / "queries"),
    target_db: Path = typer.Option(GOPC_SEARCH_DIR / "target_db" / "gopc_db"),
    results_dir: Path = typer.Option(GOPC_SEARCH_DIR / "results"),
    tmp_dir: Path = typer.Option(GOPC_SEARCH_DIR / "tmp"),
    sensitivity: float = typer.Option(7.0, "--sensitivity", "-s"),
    evalue: float = typer.Option(1e-5, "--evalue", "-e"),
    coverage: float = typer.Option(0.0, "--coverage", "-c", help="mmseqs -c; deliberately 0.0 here -- filter on qcov/tcov afterward, don't lose fragmented hits now"),
    max_seqs: int = typer.Option(10000, help="mmseqs --max-seqs -- see query_hit_cap_warning in the summary if this was too low"),
    cap_warning_fraction: float = typer.Option(0.95, help="flag a query if its hit count reaches this fraction of max-seqs"),
    mmseqs_bin: str = typer.Option("mmseqs"),
    threads: int = typer.Option(None),
    split_memory_limit: str = typer.Option(
        None,
        help="mmseqs --split-memory-limit, e.g. '50G' -- strongly recommended for a CPU (non --gpu) run under "
        "SLURM/any cgroup-limited scheduler, set comfortably below your job's --mem. Without it mmseqs sizes its "
        "prefilter split against the NODE's total memory, not your job's actual allocation, which can OOM-kill "
        "the prefilter step ('Error: Prefilter died') even though the query database itself loaded fine.",
    ),
    gpu: bool = typer.Option(
        False,
        help="mmseqs --gpu 1 -- needs a GPU-enabled mmseqs binary (mmseqs-linux-gpu, see cluster/install_mmseqs2.sh), "
        "a target database built with `gopc-build-db --gpu-compatible`, and an actual GPU available to the process. "
        "cluster/run_gopc_search.sbatch sets this up on the gpu-a100 partition.",
    ),
):
    """Sensitive MMseqs2 search of one/several/all family query FASTAs
    against the GOPC target database. Per family, writes <family>_hits.tsv
    (every alignment), <family>_unique_targets.tsv (one row per unique
    GOPC hit, best alignment retained), and <family>_summary.tsv. Does NOT
    classify hits as true PHA proteins -- alignment statistics only."""
    if family == "all":
        family_ids = sorted(p.stem for p in queries_dir.glob("*.faa"))
        if not family_ids:
            console.print(f"[red]No query FASTAs found in {queries_dir} -- run export-queries first.[/red]")
            raise typer.Exit(code=2)
    else:
        family_ids = [f.strip() for f in family.split(",") if f.strip()]

    for family_id in family_ids:
        query_fasta = queries_dir / f"{family_id}.faa"
        if not query_fasta.exists():
            console.print(f"[yellow]{query_fasta} not found -- skipping {family_id}.[/yellow]")
            continue
        console.print(f"[bold]{family_id}[/bold]: searching against GOPC (max_seqs={max_seqs})...")
        try:
            summary = gopc_search_pipeline.run_family_search(
                family_id, query_fasta, target_db, results_dir, tmp_dir,
                sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
                cap_warning_fraction=cap_warning_fraction, mmseqs_bin=mmseqs_bin, threads=threads,
                split_memory_limit=split_memory_limit, gpu=gpu,
            )
        except gopc_search_pipeline.MMseqsNotFoundError as exc:
            console.print(f"[red]{exc}[/red]")
            raise typer.Exit(code=2)

        console.print(
            f"  {summary.n_reference_queries} queries -> {summary.n_alignments} alignments -> "
            f"{summary.n_unique_gopc_targets} unique GOPC targets"
        )
        if summary.query_hit_cap_warning:
            console.print(
                f"[yellow]WARNING: query_hit_cap_warning=TRUE for {family_id} -- "
                f"{len(summary.queries_at_cap)} quer{'y' if len(summary.queries_at_cap) == 1 else 'ies'} "
                f"reached >={cap_warning_fraction:.0%} of --max-seqs {max_seqs}. "
                f"Sensitivity may be capped -- consider rerunning this family with a higher --max-seqs "
                f"(e.g. 50000/100000). See {family_id}_summary.tsv's queries_at_cap for which ones.[/yellow]"
            )


@app.command("gopc-search-batch")
def gopc_search_batch_cmd(
    family: str = typer.Option(..., help="single family_id -- run this once per family, not comma-separated/'all'"),
    n_batches: int = typer.Option(..., "--n-batches", help="total number of batches this family's queries are split into"),
    batch_index: int = typer.Option(..., "--batch-index", help="which batch this invocation runs (0-based, < n-batches)"),
    queries_dir: Path = typer.Option(GOPC_SEARCH_DIR / "queries"),
    target_db: Path = typer.Option(GOPC_SEARCH_DIR / "target_db" / "gopc_db"),
    results_dir: Path = typer.Option(GOPC_SEARCH_DIR / "results"),
    tmp_dir: Path = typer.Option(GOPC_SEARCH_DIR / "tmp"),
    sensitivity: float = typer.Option(7.0, "--sensitivity", "-s"),
    evalue: float = typer.Option(1e-5, "--evalue", "-e"),
    coverage: float = typer.Option(0.0, "--coverage", "-c"),
    max_seqs: int = typer.Option(10000, help="mmseqs --max-seqs, same value must be used for every batch of a family"),
    mmseqs_bin: str = typer.Option("mmseqs"),
    threads: int = typer.Option(None),
    split_memory_limit: str = typer.Option(None),
    gpu: bool = typer.Option(False),
):
    """Runs ONE batch of a family's queries against GOPC (see
    pipeline/gopc_search.py's run_family_search_batch docstring for why
    this exists: a single mmseqs invocation against all of GOPC restarts
    from scratch if a job is killed partway through, so a family whose
    full search doesn't fit in one job's time limit needs to be split into
    independent, individually-completable batches instead). Run once per
    batch_index in 0..n_batches-1 (a SLURM job array is the natural way --
    see cluster/run_gopc_search_batched.sbatch), then `gopc-finalize-batches`
    once every batch has a completed hits file."""
    query_fasta = queries_dir / f"{family}.faa"
    if not query_fasta.exists():
        console.print(f"[red]{query_fasta} not found -- run export-queries first.[/red]")
        raise typer.Exit(code=2)
    try:
        hits_path = gopc_search_pipeline.run_family_search_batch(
            family, query_fasta, target_db, results_dir, tmp_dir, n_batches, batch_index,
            sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
            mmseqs_bin=mmseqs_bin, threads=threads, split_memory_limit=split_memory_limit, gpu=gpu,
        )
    except gopc_search_pipeline.MMseqsNotFoundError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=2)
    console.print(f"[green]{hits_path}[/green]")


@app.command("gopc-finalize-batches")
def gopc_finalize_batches_cmd(
    family: str = typer.Option(..., help="family_id, comma-separated list, or 'all' (every family with queries/<family>.faa present)"),
    n_batches: int = typer.Option(..., "--n-batches", help="must match what gopc-search-batch was run with for this family"),
    queries_dir: Path = typer.Option(GOPC_SEARCH_DIR / "queries"),
    results_dir: Path = typer.Option(GOPC_SEARCH_DIR / "results"),
    max_seqs: int = typer.Option(10000, help="must match what gopc-search-batch was run with, used for the query_hit_cap_warning threshold"),
    cap_warning_fraction: float = typer.Option(0.95),
):
    """Concatenates every batch's hits file for a family into the same
    <family>_hits.tsv / _unique_targets.tsv / _summary.tsv a single-shot
    gopc-search run would have produced. Requires every batch
    0..n_batches-1 to have a completed hits file already (gopc-search-batch)."""
    if family == "all":
        family_ids = sorted(p.stem for p in queries_dir.glob("*.faa"))
    else:
        family_ids = [f.strip() for f in family.split(",") if f.strip()]

    for family_id in family_ids:
        query_fasta = queries_dir / f"{family_id}.faa"
        try:
            summary = gopc_search_pipeline.finalize_family_batches(
                family_id, query_fasta, results_dir, n_batches, max_seqs=max_seqs, cap_warning_fraction=cap_warning_fraction
            )
        except FileNotFoundError as exc:
            console.print(f"[red]{exc}[/red]")
            raise typer.Exit(code=2)
        console.print(
            f"[bold]{family_id}[/bold]: {summary.n_reference_queries} queries -> {summary.n_alignments} alignments -> "
            f"{summary.n_unique_gopc_targets} unique GOPC targets"
        )
        if summary.query_hit_cap_warning:
            console.print(
                f"[yellow]WARNING: query_hit_cap_warning=TRUE for {family_id} -- "
                f"{len(summary.queries_at_cap)} quer{'y' if len(summary.queries_at_cap) == 1 else 'ies'} "
                f"reached >={cap_warning_fraction:.0%} of --max-seqs {max_seqs}.[/yellow]"
            )


@app.command("gopc-combine")
def gopc_combine_cmd(results_dir: Path = typer.Option(GOPC_SEARCH_DIR / "results")):
    """Concatenates every family's *_unique_targets.tsv / *_summary.tsv
    into all_families_unique_targets.tsv / all_families_summary.tsv.
    Deliberately does NOT resolve a GOPC target hit by multiple families --
    both rows are kept."""
    n_targets, n_families = gopc_search_pipeline.combine_all_families(results_dir)
    console.print(f"[green]{results_dir / 'all_families_unique_targets.tsv'}: {n_targets} rows[/green]")
    console.print(f"[green]{results_dir / 'all_families_summary.tsv'}: {n_families} families[/green]")


@app.command("status")
def status_cmd():
    """Summary counts: proteins, family assignments by tier, phenotypes,
    retrieval runs."""
    with session_scope() as session:
        n_proteins = session.scalar(select(func.count()).select_from(Protein))
        n_evidence = session.scalar(select(func.count()).select_from(SourceEvidence))
        n_phenotypes = session.scalar(select(func.count()).select_from(Phenotype))
        n_runs = session.scalar(select(func.count()).select_from(RetrievalRun))

        console.print(f"proteins={n_proteins}  source_evidence={n_evidence}  phenotypes={n_phenotypes}  retrieval_runs={n_runs}")

        table = Table(title="family_assignment by family / evidence_tier")
        table.add_column("family_id")
        table.add_column("evidence_tier")
        table.add_column("count", justify="right")
        table.add_column("needs_manual_review", justify="right")
        rows = session.execute(
            select(
                FamilyAssignment.family_id,
                FamilyAssignment.evidence_tier,
                func.count(),
                func.sum(FamilyAssignment.needs_manual_review),
            ).group_by(FamilyAssignment.family_id, FamilyAssignment.evidence_tier)
        ).all()
        for family_id, tier, count, review_count in sorted(rows):
            table.add_row(family_id, tier, str(count), str(review_count or 0))
        console.print(table)


if __name__ == "__main__":
    app()
