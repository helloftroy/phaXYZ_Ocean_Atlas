#!/usr/bin/env python3
"""Export the 5 CBB-positive Thioglobus/SUP05 genomes for reannotation."""

from __future__ import annotations

import csv
import gzip
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PRESENCE = ROOT / "figures" / "thioglobus_sup05_kegg_marker_presence.tsv"
OUT_BASE = ROOT / "thioglobus5_reannotation"
PROTEIN_DIR = OUT_BASE / "proteins"
GENOME_DIR = OUT_BASE / "genomes"
RAW_DIR = OUT_BASE / "raw_downloads"
MANIFEST = OUT_BASE / "thioglobus5_manifest.tsv"
COMBINED_PROTEINS = OUT_BASE / "thioglobus5_combined.faa"


def api_url(genome: str, file_type: str) -> str:
    return (
        "https://motus-api.microbiomics.io/v1/genomes/"
        f"{urllib.parse.quote(genome)}"
        f"/download?file_type={urllib.parse.quote(file_type)}"
    )


def fetch(url: str, out: Path) -> str:
    if out.exists() and out.stat().st_size > 0:
        return "cached"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "pha-ocean-atlas-thioglobus5/1.0"})
    context = ssl._create_unverified_context()
    try:
        with urllib.request.urlopen(req, timeout=120, context=context) as response:
            tmp.write_bytes(response.read())
    except urllib.error.HTTPError as exc:
        return f"http_{exc.code}"
    except urllib.error.URLError as exc:
        return f"url_error:{exc.reason}"
    if tmp.stat().st_size == 0:
        tmp.unlink(missing_ok=True)
        return "empty"
    tmp.replace(out)
    time.sleep(0.15)
    return "downloaded"


def gunzip_to_text(src: Path, dst: Path) -> int:
    n_records = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rt", errors="replace") as inp, dst.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                n_records += 1
            out.write(line)
    return n_records


def main() -> int:
    with PRESENCE.open() as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["has_cbb_rbcLS_prk"] == "yes"
        ]
    if len(rows) != 5:
        raise SystemExit(f"Expected 5 CBB-positive Thioglobus/SUP05 genomes; found {len(rows)}")

    manifest_rows = []
    for row in rows:
        genome = row["genome"]
        raw_faa = RAW_DIR / f"{genome}.genes.faa.gz"
        raw_genome = RAW_DIR / f"{genome}.fa.gz"
        faa_status = fetch(api_url(genome, "gene_faa"), raw_faa)
        genome_status = fetch(api_url(genome, "genome"), raw_genome)

        faa = PROTEIN_DIR / f"{genome}.faa"
        genome_fasta = GENOME_DIR / f"{genome}.fasta"
        n_proteins = gunzip_to_text(raw_faa, faa) if raw_faa.exists() else 0
        n_contigs = gunzip_to_text(raw_genome, genome_fasta) if raw_genome.exists() else 0

        manifest_rows.append(
            {
                "genome": genome,
                "gtdb_species": row["gtdb_species"],
                "omdb_completeness_pct": row["omdb_completeness_pct"],
                "omdb_contamination_pct": row["omdb_contamination_pct"],
                "n_proteins": n_proteins,
                "n_contigs": n_contigs,
                "protein_fasta": str(faa.relative_to(ROOT)),
                "genome_fasta": str(genome_fasta.relative_to(ROOT)),
                "protein_download_status": faa_status,
                "genome_download_status": genome_status,
            }
        )

    with COMBINED_PROTEINS.open("w") as out:
        for row in manifest_rows:
            out.write((ROOT / row["protein_fasta"]).read_text())

    fields = [
        "genome",
        "gtdb_species",
        "omdb_completeness_pct",
        "omdb_contamination_pct",
        "n_proteins",
        "n_contigs",
        "protein_fasta",
        "genome_fasta",
        "protein_download_status",
        "genome_download_status",
    ]
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"Wrote {MANIFEST.relative_to(ROOT)}")
    print(f"Wrote {COMBINED_PROTEINS.relative_to(ROOT)}")
    print(f"Protein FASTAs: {PROTEIN_DIR.relative_to(ROOT)}")
    print(f"Genome FASTAs: {GENOME_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
