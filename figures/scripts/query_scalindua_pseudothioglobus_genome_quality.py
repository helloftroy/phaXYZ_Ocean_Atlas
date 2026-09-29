#!/usr/bin/env python3
"""Genome-quality follow-up for the Scalindua phaC finding (§12.1/12.5):
BEMA21-1_SAMN15000289_MAG_00000040 (Scalindua, the only phaC-positive
Scalindua genome) and GLAS15-1_SAMN02905564_MAG_00000008 (Pseudothioglobus)
share a byte-identical NR100 phaC sequence (OMDBv2.0_AA_G_NR100_000013537294,
the phaC_cluster0.7 representative for a 38-genome, mostly-Pseudothioglobus
cluster). completeness/contamination for both are already cached in
phaC_unique_targets_with_metadata.tsv; this re-fetches the same two
numbers live from the OMDB genome-browser page (same window._MB_GENOME_DATA
JSON block query_thioglobus_sup05_kegg_markers.py already parses) purely
as an independent cross-check, plus pulls the extra anvio/checkm detail
and assembly stats (genome_size, n50, no_scaffolds) that
phaC_unique_targets_with_metadata.tsv does not carry.
"""
from __future__ import annotations

import csv
import json
import re
import ssl
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "scalindua_pseudothioglobus_genome_quality.tsv"
GENOME_JSON_RE = re.compile(r"window\._MB_GENOME_DATA = JSON\.parse\(`(.+?)`\);", re.S)

GENOMES = [
    "BEMA21-1_SAMN15000289_MAG_00000040",
    "GLAS15-1_SAMN02905564_MAG_00000008",
]

FIELDS = [
    "genome", "genus", "species", "is_mag", "completeness", "contamination",
    "checkm_completeness", "checkm_contamination", "anvio_completion", "anvio_contamination",
    "genome_size", "n50", "no_scaffolds", "bioproject", "biosample", "publication",
]


def fetch_genome_metadata(genome: str) -> dict[str, object]:
    url = f"https://omdb.microbiomics.io/repository/ocean/genome/{urllib.parse.quote(genome)}"
    req = urllib.request.Request(url, headers={"User-Agent": "pha-ocean-atlas-genome-quality-query/1.0"})
    context = ssl._create_unverified_context()
    with urllib.request.urlopen(req, timeout=60, context=context) as resp:
        text = resp.read().decode(errors="replace")
    match = GENOME_JSON_RE.search(text)
    if not match:
        return {}
    return json.loads(match.group(1))


def main() -> int:
    rows = []
    for genome in GENOMES:
        data = fetch_genome_metadata(genome)
        rows.append({
            "genome": genome,
            "genus": data.get("genus"),
            "species": data.get("species"),
            "is_mag": data.get("is_mag"),
            "completeness": data.get("completeness"),
            "contamination": data.get("contamination"),
            "checkm_completeness": data.get("checkm_completeness"),
            "checkm_contamination": data.get("checkm_contamination"),
            "anvio_completion": data.get("anvio_completion"),
            "anvio_contamination": data.get("anvio_contamination"),
            "genome_size": data.get("genome_size"),
            "n50": data.get("n50"),
            "no_scaffolds": data.get("no_scaffolds"),
            "bioproject": data.get("bioproject"),
            "biosample": data.get("public_sample_link", "").rsplit("/", 1)[-1] if data.get("public_sample_link") else "",
            "publication": data.get("publication"),
        })
        print(genome, "->", data.get("completeness"), "% complete,", data.get("contamination"), "% contamination")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
