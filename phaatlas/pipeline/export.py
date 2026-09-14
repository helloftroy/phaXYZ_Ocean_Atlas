"""Generates PHA_reference/exports/pha_reference_master.csv and
pha_reference.faa FROM the SQLite database -- never the other way around.
Both are fully regenerable from phaatlas.sqlite at any time; neither
should ever be hand-edited (see README).
"""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

from phaatlas.db.session import connect_readonly

MASTER_COLUMNS = [
    "protein_id",
    "pha_family",
    "pha_subfamily",
    "uniprot_accession",
    "secondary_accessions",
    "reviewed",
    "protein_name",
    "alternative_protein_names",
    "gene_name",
    "gene_synonyms",
    "organism",
    "strain",
    "ncbi_taxid",
    "taxonomic_lineage",
    "sequence_length",
    "sequence_sha256",
    "sequence",
    "protein_existence",
    "annotation_score",
    "ec_numbers",
    "rhea_ids",
    "function_text",
    "catalytic_activity_text",
    "pathway_text",
    "brenda_ec_numbers",
    "brenda_record_ids",
    "brenda_sources",
    "retrieval_reason",
    "evidence_tier",
    "cluster95_id",
    "cluster95_representative",
    "cluster95_size",
    "literature_pmids",
    "literature_dois",
    "demonstrated_monomers",
    "demonstrated_polymers",
    "polymer_composition_text",
    "substrate_specificity_text",
    "phenotype_evidence_level",
    "phenotype_references",
    "temperature_optimum",
    "temperature_range",
    "temperature_stability",
    "ph_optimum",
    "ph_range",
    "specific_activity",
    "source_releases",
    "retrieved_at",
]

# Columns where duplicate pipe-separated tokens can legitimately accumulate
# across multiple evidence rows feeding the same view row (e.g. the same
# PMID cited by two different BRENDA substrate entries) -- deduped here
# (order-preserving) so the CSV itself stays clean, without complicating the
# SQL view with per-token DISTINCT logic SQLite can't express directly.
DEDUP_PIPE_COLUMNS = {
    "secondary_accessions",
    "alternative_protein_names",
    "gene_synonyms",
    "ec_numbers",
    "rhea_ids",
    "brenda_ec_numbers",
    "brenda_record_ids",
    "brenda_sources",
    "retrieval_reason",
    "literature_pmids",
    "literature_dois",
    "demonstrated_monomers",
    "demonstrated_polymers",
    "phenotype_evidence_level",
    "phenotype_references",
    "source_releases",
}


def _dedup_pipe(value: str | None) -> str | None:
    if not value:
        return value
    seen: list[str] = []
    for token in value.split("|"):
        if token and token not in seen:
            seen.append(token)
    return "|".join(seen)


def export_master_csv(db_path: Path, output_csv_path: Path) -> int:
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
    conn.row_factory = sqlite3.Row
    try:
        cursor = conn.execute(f"SELECT {', '.join(MASTER_COLUMNS)} FROM protein_master_export")
        rows = cursor.fetchall()
    finally:
        conn.close()

    with open(output_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MASTER_COLUMNS)
        writer.writeheader()
        for row in rows:
            record = dict(row)
            for col in DEDUP_PIPE_COLUMNS:
                record[col] = _dedup_pipe(record.get(col))
            writer.writerow(record)

    return len(rows)


def export_fasta(db_path: Path, output_faa_path: Path) -> int:
    output_faa_path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
    conn.row_factory = sqlite3.Row
    try:
        cursor = conn.execute(
            """
            SELECT protein_id, uniprot_accession, gene_name, organism, sequence
            FROM protein
            WHERE sequence IS NOT NULL AND sequence != ''
            ORDER BY protein_id
            """
        )
        rows = cursor.fetchall()
    finally:
        conn.close()

    count = 0
    with open(output_faa_path, "w") as f:
        for row in rows:
            header_bits = [row["protein_id"]]
            if row["gene_name"]:
                header_bits.append(f"gene={row['gene_name']}")
            if row["organism"]:
                header_bits.append(f"organism={row['organism']}")
            f.write(">" + " ".join(header_bits) + "\n")
            seq = row["sequence"]
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
            count += 1

    return count


def _write_fasta_records(rows, output_path: Path) -> int:
    count = 0
    with open(output_path, "w") as f:
        for row in rows:
            header_bits = [row["protein_id"], f"family={row['pha_family']}"]
            if row["gene_name"]:
                header_bits.append(f"gene={row['gene_name']}")
            if row["organism"]:
                header_bits.append(f"organism={row['organism']}")
            if "cluster95_size" in row.keys() and row["cluster95_size"] is not None:
                header_bits.append(f"cluster95_size={row['cluster95_size']}")
            f.write(">" + " ".join(header_bits) + "\n")
            seq = row["sequence"]
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")
            count += 1
    return count


def export_fasta_all(db_path: Path, output_faa_path: Path) -> int:
    """pha_reference_all.faa -- every (protein, family) row with a sequence,
    straight off protein_master_export, headers carrying family/gene/
    organism. Unlike export_fasta() (the bare `protein` table, one entry
    per accession), this is family-scoped like the rest of the cluster95
    bundle -- a protein assigned to two families appears twice, once per
    family, matching the master CSV's own row granularity.
    """
    output_faa_path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            """
            SELECT protein_id, pha_family, gene_name, organism, sequence
            FROM protein_master_export
            WHERE sequence IS NOT NULL AND sequence != '' AND pha_family IS NOT NULL
            ORDER BY pha_family, protein_id
            """
        ).fetchall()
    finally:
        conn.close()
    return _write_fasta_records(rows, output_faa_path)


def export_fasta_nr95(db_path: Path, output_faa_path: Path) -> int:
    """pha_reference_nr95.faa -- one sequence per (family, cluster95_id):
    exactly the rows whose protein_id IS that cluster's cluster95_representative.
    Requires `pha-reference cluster95` to have populated cluster95_* first;
    a family/protein that hasn't been clustered yet (or has no sequence)
    simply doesn't appear here.
    """
    output_faa_path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            """
            SELECT protein_id, pha_family, gene_name, organism, sequence, cluster95_size
            FROM protein_master_export
            WHERE sequence IS NOT NULL AND sequence != ''
              AND cluster95_representative IS NOT NULL
              AND protein_id = cluster95_representative
            ORDER BY pha_family, protein_id
            """
        ).fetchall()
    finally:
        conn.close()
    return _write_fasta_records(rows, output_faa_path)


def export_cluster95_tsv(db_path: Path, output_tsv_path: Path) -> int:
    """pha_reference_cluster95.tsv -- flat membership table, one row per
    (protein, family) that has been clustered, independent of SQLite for
    downstream tools that would rather not open the database directly.
    """
    output_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect_readonly(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            """
            SELECT protein_id, pha_family, cluster95_id, cluster95_representative, cluster95_size
            FROM protein_master_export
            WHERE cluster95_id IS NOT NULL
            ORDER BY pha_family, cluster95_id, protein_id
            """
        ).fetchall()
    finally:
        conn.close()

    with open(output_tsv_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id", "pha_family", "cluster95_id", "cluster95_representative", "cluster95_size"])
        for row in rows:
            writer.writerow(list(row))
    return len(rows)
