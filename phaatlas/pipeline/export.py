"""Generates PHA_reference/exports/pha_reference_master.csv and
pha_reference.faa FROM the SQLite database -- never the other way around.
Both are fully regenerable from phaatlas.sqlite at any time; neither
should ever be hand-edited (see README).
"""

from __future__ import annotations

import csv
import sqlite3
from pathlib import Path

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
    conn = sqlite3.connect(str(db_path))
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
    conn = sqlite3.connect(str(db_path))
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
