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
from phaatlas.pipeline import depth_heatmap as depth_heatmap_pipeline
from phaatlas.pipeline import gopc_search as gopc_search_pipeline
from phaatlas.pipeline import ncbi_depth as ncbi_depth_pipeline
from phaatlas.pipeline import omdb_metadata as omdb_metadata_pipeline
from phaatlas.pipeline import pathway_architecture as pathway_architecture_pipeline
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
OMDB_SEARCH_DIR = REPO_ROOT / "PHA_bioprospecting" / "omdb_search"
# Queries are the same NR95 reference FASTAs regardless of which catalog
# they're searched against -- shared/reused rather than duplicated. Kept
# under gopc_search/ (not its own directory) because that path is already
# hardcoded into the currently-running GOPC batch jobs' sbatch defaults;
# moving it would break those mid-flight for a purely cosmetic gain.
SHARED_QUERIES_DIR = GOPC_SEARCH_DIR / "queries"
OMDB_DATABASES_DIR = REPO_ROOT / "PHA_bioprospecting" / "databases" / "OMDBv2"


@app.command("export-queries")
def export_queries_cmd(
    queries_dir: Path = typer.Option(SHARED_QUERIES_DIR, help="output directory, one <family>.faa per family -- shared between GOPC and OMDB searches"),
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
    db_load_mode: int = typer.Option(
        None,
        help="mmseqs --db-load-mode (0 auto / 1 fread / 2 mmap / 3 mmap+touch). Confirmed live: convertalis (the "
        "final output-writing step) has no dedicated memory cap the way --split-memory-limit covers the prefilter, "
        "and can still OOM-kill a large family even with a generous --mem. Explicit mmap (2) is the untested next "
        "lever to try if that happens -- cgroups can reclaim clean mmap'd pages under pressure before OOM-killing.",
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
                split_memory_limit=split_memory_limit, gpu=gpu, db_load_mode=db_load_mode,
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
    db_load_mode: int = typer.Option(None, help="mmseqs --db-load-mode -- see gopc-search --help for the full explanation"),
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
            db_load_mode=db_load_mode,
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


# --- OMDB commands: identical pipeline to the gopc-* commands above (every
# function in pipeline/gopc_search.py takes target_db_path as a plain
# parameter, nothing GOPC-specific baked in) -- just pointed at
# OMDB_SEARCH_DIR instead of GOPC_SEARCH_DIR, and sharing the same
# queries/ directory (see SHARED_QUERIES_DIR above). Kept as separate
# omdb-* commands rather than a --catalog flag on the gopc-* ones so the
# two searches' results/target_db/tmp never collide and neither command
# set risks the other mid-flight.


@app.command("omdb-build-db")
def omdb_build_db_cmd(
    omdb_faa: Path = typer.Argument(..., help="OMDBv2.0_AA_G_NR100.faa.gz, .gz is fine -- mmseqs reads gzip directly"),
    target_db: Path = typer.Option(OMDB_SEARCH_DIR / "target_db" / "omdb_db", help="mmseqs target database path (prefix)"),
    mmseqs_bin: str = typer.Option("mmseqs"),
    threads: int = typer.Option(None),
    shuffle: bool = typer.Option(False, help="see gopc-build-db --shuffle -- same reasoning applies"),
    gpu_compatible: bool = typer.Option(False, help="see gopc-build-db --gpu-compatible -- required before omdb-search --gpu"),
):
    """`mmseqs createdb` on OMDBv2.0_AA_G_NR100 -- run once. ~10x fewer
    sequences than GOPC (249.5M vs. ~2.46B), so this and the search itself
    should need noticeably less time/memory -- but that's an expectation,
    not something we've confirmed live yet."""
    try:
        gopc_search_pipeline.build_gopc_target_db(
            omdb_faa, target_db, mmseqs_bin=mmseqs_bin, threads=threads, shuffle=shuffle, gpu_compatible=gpu_compatible
        )
    except gopc_search_pipeline.MMseqsNotFoundError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=2)
    console.print(f"[green]OMDB target database ready at {target_db}[/green]")


@app.command("omdb-search")
def omdb_search_cmd(
    family: str = typer.Option("all", help="family_id, comma-separated list, or 'all'"),
    queries_dir: Path = typer.Option(SHARED_QUERIES_DIR),
    target_db: Path = typer.Option(OMDB_SEARCH_DIR / "target_db" / "omdb_db"),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results"),
    tmp_dir: Path = typer.Option(OMDB_SEARCH_DIR / "tmp"),
    sensitivity: float = typer.Option(7.0, "--sensitivity", "-s"),
    evalue: float = typer.Option(1e-5, "--evalue", "-e"),
    coverage: float = typer.Option(0.0, "--coverage", "-c"),
    max_seqs: int = typer.Option(10000),
    cap_warning_fraction: float = typer.Option(0.95),
    mmseqs_bin: str = typer.Option("mmseqs"),
    threads: int = typer.Option(None),
    split_memory_limit: str = typer.Option(None),
    gpu: bool = typer.Option(False),
    db_load_mode: int = typer.Option(None, help="mmseqs --db-load-mode -- see gopc-search --help for the full explanation"),
):
    """Same as gopc-search, against OMDBv2.0_AA_G_NR100 instead of GOPC.
    Does NOT classify hits as true PHA proteins -- alignment statistics only."""
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
        console.print(f"[bold]{family_id}[/bold]: searching against OMDB (max_seqs={max_seqs})...")
        try:
            summary = gopc_search_pipeline.run_family_search(
                family_id, query_fasta, target_db, results_dir, tmp_dir,
                sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
                cap_warning_fraction=cap_warning_fraction, mmseqs_bin=mmseqs_bin, threads=threads,
                split_memory_limit=split_memory_limit, gpu=gpu, db_load_mode=db_load_mode,
            )
        except gopc_search_pipeline.MMseqsNotFoundError as exc:
            console.print(f"[red]{exc}[/red]")
            raise typer.Exit(code=2)

        console.print(
            f"  {summary.n_reference_queries} queries -> {summary.n_alignments} alignments -> "
            f"{summary.n_unique_gopc_targets} unique OMDB targets"
        )
        if summary.query_hit_cap_warning:
            console.print(
                f"[yellow]WARNING: query_hit_cap_warning=TRUE for {family_id} -- "
                f"{len(summary.queries_at_cap)} quer{'y' if len(summary.queries_at_cap) == 1 else 'ies'} "
                f"reached >={cap_warning_fraction:.0%} of --max-seqs {max_seqs}.[/yellow]"
            )


@app.command("omdb-search-batch")
def omdb_search_batch_cmd(
    family: str = typer.Option(..., help="single family_id"),
    n_batches: int = typer.Option(..., "--n-batches"),
    batch_index: int = typer.Option(..., "--batch-index"),
    queries_dir: Path = typer.Option(SHARED_QUERIES_DIR),
    target_db: Path = typer.Option(OMDB_SEARCH_DIR / "target_db" / "omdb_db"),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results"),
    tmp_dir: Path = typer.Option(OMDB_SEARCH_DIR / "tmp"),
    sensitivity: float = typer.Option(7.0, "--sensitivity", "-s"),
    evalue: float = typer.Option(1e-5, "--evalue", "-e"),
    coverage: float = typer.Option(0.0, "--coverage", "-c"),
    max_seqs: int = typer.Option(10000),
    mmseqs_bin: str = typer.Option("mmseqs"),
    threads: int = typer.Option(None),
    split_memory_limit: str = typer.Option(None),
    gpu: bool = typer.Option(False),
    db_load_mode: int = typer.Option(None, help="mmseqs --db-load-mode -- see gopc-search --help for the full explanation"),
):
    """Same as gopc-search-batch, against OMDBv2.0_AA_G_NR100 instead of
    GOPC -- included for parity/in case OMDB turns out to need it too, even
    though it's expected to be much less likely given OMDB's smaller size."""
    query_fasta = queries_dir / f"{family}.faa"
    if not query_fasta.exists():
        console.print(f"[red]{query_fasta} not found -- run export-queries first.[/red]")
        raise typer.Exit(code=2)
    try:
        hits_path = gopc_search_pipeline.run_family_search_batch(
            family, query_fasta, target_db, results_dir, tmp_dir, n_batches, batch_index,
            sensitivity=sensitivity, evalue=evalue, coverage=coverage, max_seqs=max_seqs,
            mmseqs_bin=mmseqs_bin, threads=threads, split_memory_limit=split_memory_limit, gpu=gpu,
            db_load_mode=db_load_mode,
        )
    except gopc_search_pipeline.MMseqsNotFoundError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=2)
    console.print(f"[green]{hits_path}[/green]")


@app.command("omdb-finalize-batches")
def omdb_finalize_batches_cmd(
    family: str = typer.Option(..., help="family_id, comma-separated list, or 'all'"),
    n_batches: int = typer.Option(..., "--n-batches"),
    queries_dir: Path = typer.Option(SHARED_QUERIES_DIR),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results"),
    max_seqs: int = typer.Option(10000),
    cap_warning_fraction: float = typer.Option(0.95),
):
    """Same as gopc-finalize-batches, for OMDB."""
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
            f"{summary.n_unique_gopc_targets} unique OMDB targets"
        )
        if summary.query_hit_cap_warning:
            console.print(f"[yellow]WARNING: query_hit_cap_warning=TRUE for {family_id}[/yellow]")


@app.command("omdb-combine")
def omdb_combine_cmd(results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results")):
    """Same as gopc-combine, for OMDB."""
    n_targets, n_families = gopc_search_pipeline.combine_all_families(results_dir)
    console.print(f"[green]{results_dir / 'all_families_unique_targets.tsv'}: {n_targets} rows[/green]")
    console.print(f"[green]{results_dir / 'all_families_summary.tsv'}: {n_families} families[/green]")


@app.command("omdb-enrich-metadata")
def omdb_enrich_metadata_cmd(
    family: str = typer.Option(..., help="single family_id or comma-separated list -- NOT 'all', deliberately: "
                                          "start with a small family (e.g. phaQ, ~5K rows) before a large one "
                                          "(phaB/phaA are 500K+), since cost scales with distinct genomes found"),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results"),
    cluster_tsv: Path = typer.Option(OMDB_DATABASES_DIR / "OMDBv2.0_AA_G_NR100.cluster.tsv.gz",
                                      help="OMDBv2.0_AA_G_NR100.cluster.tsv.gz -- download via "
                                           "./PHA_bioprospecting/scripts/download_omdb.sh nr100-clusters"),
    max_genomes_per_target: int = typer.Option(5, help="cap on how many distinct genomes to report per NR100 "
                                                         "cluster -- n_genomes_in_cluster_total is never capped"),
    genome_batch_size: int = typer.Option(100),
    sample_batch_size: int = typer.Option(100),
    api_base: str = typer.Option(omdb_metadata_pipeline.OMDB_API_BASE),
    sleep_seconds: float = typer.Option(0.3, help="politeness delay between batched API requests"),
):
    """Joins a family's <family>_unique_targets.tsv to OMDB's own genome
    (GTDB taxonomy) and sample (lat/lon, ecosystem) metadata, writing
    <family>_unique_targets_with_metadata.tsv. Needs internet (OMDB's
    public API) and cluster_tsv downloaded locally -- neither GPU nor
    mmseqs required. See pipeline/omdb_metadata.py's module docstring for
    how target_id -> genome -> sample -> metadata is resolved."""
    if not cluster_tsv.exists():
        console.print(f"[red]{cluster_tsv} not found -- download it first: "
                       f"./PHA_bioprospecting/scripts/download_omdb.sh nr100-clusters[/red]")
        raise typer.Exit(code=2)

    family_ids = [f.strip() for f in family.split(",") if f.strip()]
    for family_id in family_ids:
        unique_targets_path = results_dir / f"{family_id}_unique_targets.tsv"
        if not unique_targets_path.exists():
            console.print(f"[yellow]{unique_targets_path} not found -- skipping {family_id}.[/yellow]")
            continue
        out_path = results_dir / f"{family_id}_unique_targets_with_metadata.tsv"
        console.print(f"[bold]{family_id}[/bold]: resolving target_id -> genome via {cluster_tsv}...")
        summary = omdb_metadata_pipeline.enrich_unique_targets(
            unique_targets_path, cluster_tsv, out_path,
            max_genomes_per_target=max_genomes_per_target,
            genome_batch_size=genome_batch_size, sample_batch_size=sample_batch_size,
            api_base=api_base, sleep_seconds=sleep_seconds,
        )
        console.print(
            f"  {summary.n_target_ids} target_ids -> {summary.n_targets_matched_in_cluster_file} matched in cluster file\n"
            f"  {summary.n_distinct_genomes_needed} distinct genomes needed -> {summary.n_genomes_resolved} resolved via OMDB API\n"
            f"  {summary.n_distinct_samples_needed} distinct samples needed -> {summary.n_samples_resolved} resolved via OMDB API\n"
            f"  -> [green]{out_path}[/green]: {summary.n_output_rows} rows"
        )
        if summary.n_targets_matched_in_cluster_file < summary.n_target_ids:
            missing = summary.n_target_ids - summary.n_targets_matched_in_cluster_file
            console.print(f"[yellow]WARNING: {missing} target_id(s) not found in {cluster_tsv} -- "
                           f"unexpected unless target_id came from a different NR100 release.[/yellow]")


@app.command("omdb-enrich-depth")
def omdb_enrich_depth_cmd(
    family: str = typer.Option(..., help="single family_id or comma-separated list"),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results",
                                      help="where <family>_unique_targets_with_metadata.tsv already lives "
                                           "(from omdb-enrich-metadata)"),
    batch_size: int = typer.Option(100, help="BioSample accessions per NCBI efetch request"),
    sleep_seconds: float = typer.Option(0.35, help="politeness delay between NCBI requests"),
):
    """Adds depth_raw/depth_m/depth_zone columns to an already-built
    <family>_unique_targets_with_metadata.tsv, by fetching each row's
    sample's real NCBI BioSample record (OMDB's own API has no depth
    field at all -- see pipeline/ncbi_depth.py's module docstring).
    Writes <family>_unique_targets_with_metadata_depth.tsv. Coverage is
    inherently partial -- not every sample has a real NCBI accession, and
    not every accession's record reports depth -- printed counts make
    that visible rather than silent."""
    family_ids = [f.strip() for f in family.split(",") if f.strip()]
    for family_id in family_ids:
        metadata_path = results_dir / f"{family_id}_unique_targets_with_metadata.tsv"
        if not metadata_path.exists():
            console.print(f"[yellow]{metadata_path} not found -- run omdb-enrich-metadata first. Skipping {family_id}.[/yellow]")
            continue
        out_path = results_dir / f"{family_id}_unique_targets_with_metadata_depth.tsv"
        console.print(f"[bold]{family_id}[/bold]: fetching BioSample depth records from NCBI...")
        stats = ncbi_depth_pipeline.enrich_with_depth(metadata_path, out_path, batch_size=batch_size, sleep_seconds=sleep_seconds)
        console.print(
            f"  {stats['n_distinct_biosamples']} distinct real BioSample accessions -> "
            f"{stats['n_biosamples_found_in_ncbi']} found in NCBI -> "
            f"{stats['n_rows_with_parsed_depth']}/{stats['n_input_rows']} rows got a parsed depth\n"
            f"  -> [green]{out_path}[/green]"
        )


@app.command("depth-heatmap")
def depth_heatmap_cmd(
    family: str = typer.Option(..., help="single family_id"),
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results",
                                      help="where <family>_unique_targets_with_metadata_depth.tsv lives "
                                           "(from omdb-enrich-depth)"),
    clade_column: str = typer.Option("best_query", help="which column groups rows into a 'clade' -- "
                                                          "best_query (closest NR95 reference protein) by default"),
    top_n: int = typer.Option(15, help="how many top clades to print to the console (the written TSVs are not truncated)"),
):
    """Builds a clade x depth-zone abundance matrix for one family, e.g.
    "is one PhaC clade almost exclusively deep-ocean while another
    dominates the photic zone". Writes both a long-format TSV (ready for
    pandas/seaborn) and a wide clade x depth-zone matrix TSV, ready to
    paste straight into a heatmap. See pipeline/depth_heatmap.py's module
    docstring for exactly what "clade" means here (a proxy, not a real
    phylogeny)."""
    metadata_depth_path = results_dir / f"{family}_unique_targets_with_metadata_depth.tsv"
    if not metadata_depth_path.exists():
        console.print(f"[red]{metadata_depth_path} not found -- run omdb-enrich-depth first.[/red]")
        raise typer.Exit(code=2)

    clades_ranked, genome_counts, target_id_counts = depth_heatmap_pipeline.build_clade_depth_matrix(
        metadata_depth_path, clade_column=clade_column
    )
    long_path = results_dir / f"{family}_clade_depth_long.tsv"
    wide_path = results_dir / f"{family}_clade_depth_matrix.tsv"
    n_long_rows = depth_heatmap_pipeline.write_long_format(clades_ranked, genome_counts, target_id_counts, long_path)
    depth_heatmap_pipeline.write_wide_matrix(clades_ranked, genome_counts, wide_path)

    console.print(f"[green]{long_path}[/green]: {n_long_rows} (clade, depth_zone) rows")
    console.print(f"[green]{wide_path}[/green]: {len(clades_ranked)} clades x {len(depth_heatmap_pipeline.DEPTH_ZONE_ORDER)} depth zones")

    table = Table(title=f"{family}: top clades by total genome abundance across depth zones")
    table.add_column("clade")
    for zone in depth_heatmap_pipeline.DEPTH_ZONE_ORDER:
        table.add_column(zone, justify="right")
    for clade in clades_ranked[:top_n]:
        row_counts = genome_counts.get(clade, {})
        table.add_row(clade, *(str(row_counts.get(zone, 0)) for zone in depth_heatmap_pipeline.DEPTH_ZONE_ORDER))
    console.print(table)


@app.command("pathway-architecture")
def pathway_architecture_cmd(
    results_dir: Path = typer.Option(OMDB_SEARCH_DIR / "results",
                                      help="directory to scan for *_unique_targets_with_metadata.tsv -- "
                                           "point this at wherever those files actually are, e.g. a "
                                           "scratch dir you scp'd results into"),
    out_dir: Path = typer.Option(None, help="defaults to results_dir"),
    top_n: int = typer.Option(5, help="how many top genera/species/studies to list per architecture"),
):
    """Builds a genome x PHA-family count matrix from every
    <family>_unique_targets_with_metadata.tsv found in results_dir, then
    collapses each genome's family set into a short architecture label
    (e.g. phaA+phaB+phaC -> "ABC"). Writes genome_family_matrix.tsv and
    architecture_summary.tsv, and prints the most/least common
    architectures with their top taxa. Only families that have actually
    been through omdb-enrich-metadata are reflected -- a family missing
    from results_dir reads as absent everywhere, not unknown."""
    metadata_paths = pathway_architecture_pipeline.discover_metadata_files(results_dir)
    if not metadata_paths:
        console.print(f"[red]No *_unique_targets_with_metadata.tsv found in {results_dir} -- "
                       f"run omdb-enrich-metadata first.[/red]")
        raise typer.Exit(code=2)

    families_found = sorted(p.name.removesuffix("_unique_targets_with_metadata.tsv") for p in metadata_paths)
    console.print(f"Found metadata for {len(metadata_paths)} famil{'y' if len(metadata_paths) == 1 else 'ies'}: "
                  f"{', '.join(families_found)}")
    all_families = pathway_architecture_pipeline.default_family_order()
    missing = [f for f in all_families if f not in families_found]
    if missing:
        console.print(f"[yellow]Not yet enriched (will read as absent in every architecture): {', '.join(missing)}[/yellow]")

    genomes = pathway_architecture_pipeline.load_genome_records(metadata_paths)
    out_dir = out_dir or results_dir
    matrix_path = out_dir / "genome_family_matrix.tsv"
    summary_path = out_dir / "architecture_summary.tsv"
    n_rows = pathway_architecture_pipeline.write_genome_family_matrix(genomes, all_families, matrix_path)
    stats = pathway_architecture_pipeline.summarize_architectures(genomes, all_families, top_n=top_n)
    pathway_architecture_pipeline.write_architecture_summary(stats, summary_path)

    console.print(f"[green]{matrix_path}[/green]: {n_rows} genomes")
    console.print(f"[green]{summary_path}[/green]: {len(stats)} distinct architectures")

    table = Table(title="Pathway architectures (most to least common)")
    table.add_column("architecture")
    table.add_column("n_genomes", justify="right")
    table.add_column("% of genomes", justify="right")
    table.add_column("top genus")
    table.add_column("top study")
    for s in stats:
        top_genus = f"{s.top_genera[0][0]} ({s.top_genera[0][1]})" if s.top_genera else ""
        top_study = f"{s.top_studies[0][0]} ({s.top_studies[0][1]})" if s.top_studies else ""
        table.add_row(s.architecture or "(none)", str(s.n_genomes), f"{s.pct_of_genomes:.1f}%", top_genus, top_study)
    console.print(table)


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
