from __future__ import annotations

import json
from pathlib import Path

import typer
from dotenv import load_dotenv
from rich.console import Console
from rich.table import Table
from sqlalchemy import func, select

from phaatlas.config_loader import DEFAULT_FAMILY_CONFIG_PATH, load_family_definitions
from phaatlas.db.models import FamilyAssignment, FamilyDefinitionRow, Phenotype, Protein, RetrievalRun, SourceEvidence
from phaatlas.db.session import get_db_path, init_db, session_scope
from phaatlas.pipeline import export as export_pipeline
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
