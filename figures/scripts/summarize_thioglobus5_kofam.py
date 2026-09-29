#!/usr/bin/env python3
"""Summarize KOfamScan marker calls for the 5 Thioglobus/SUP05 genomes."""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "thioglobus5_reannotation" / "thioglobus5_manifest.tsv"
KOFAM_DIR = ROOT / "thioglobus5_reannotation" / "kofam"
OUT = ROOT / "figures" / "thioglobus5_kofam_marker_summary.tsv"
HITS_OUT = ROOT / "figures" / "thioglobus5_kofam_marker_hits.tsv"

MARKERS = {
    "cbbM_rbcL": {"K01601"},
    "rbcS": {"K01602"},
    "prk": {"K00855"},
    "sqr": {"K17218"},
    "fccA": {"K17230"},
    "fccB": {"K17229"},
    "soxX": {"K17223"},
    "soxA": {"K17222"},
    "soxB": {"K17224"},
    "soxY": {"K17226"},
    "soxZ": {"K17227"},
    "soxC": {"K17225"},
    "soxD": {"K22622"},
    "dsrA": {"K11180"},
    "dsrB": {"K11181"},
    "dsrE": {"K11182"},
    "dsrF": {"K11183"},
    "dsrH": {"K11184"},
    "dsrC": {"K11179"},
    "aprA": {"K00394"},
    "aprB": {"K00395"},
    "sat": {"K00958"},
}
KO_TO_MARKERS = {ko: marker for marker, kos in MARKERS.items() for ko in kos}
SOX_CORE = {"soxX", "soxA", "soxB", "soxY", "soxZ"}
CBBCORE = {"cbbM_rbcL", "rbcS", "prk"}
KO_RE = re.compile(r"K\d{5}")


def parse_detail_tsv(path: Path) -> list[dict[str, str]]:
    hits = []
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            significant = line.startswith("*")
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if significant:
                parts[0] = parts[0].lstrip("*").strip()
            if len(parts) < 2:
                continue
            gene = parts[0].strip()
            ko_match = KO_RE.search(line)
            if not ko_match:
                continue
            ko = ko_match.group(0)
            marker = KO_TO_MARKERS.get(ko)
            if marker and significant:
                hits.append(
                    {
                        "query": gene,
                        "marker": marker,
                        "ko": ko,
                        "raw_line": line,
                    }
                )
    return hits


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def main() -> int:
    with MANIFEST.open() as handle:
        genomes = list(csv.DictReader(handle, delimiter="\t"))

    summary_rows = []
    hit_rows = []
    for genome_row in genomes:
        genome = genome_row["genome"]
        kofam_path = KOFAM_DIR / f"{genome}.kofam.detail.tsv"
        hits = parse_detail_tsv(kofam_path) if kofam_path.exists() else []
        by_marker = defaultdict(list)
        for hit in hits:
            by_marker[hit["marker"]].append(hit)
            hit_rows.append(
                {
                    "genome": genome,
                    "gtdb_species": genome_row["gtdb_species"],
                    "marker": hit["marker"],
                    "ko": hit["ko"],
                    "query": hit["query"],
                    "raw_line": hit["raw_line"],
                }
            )

        present = set(by_marker)
        out = {
            "genome": genome,
            "gtdb_species": genome_row["gtdb_species"],
            "omdb_completeness_pct": genome_row["omdb_completeness_pct"],
            "n_proteins": genome_row["n_proteins"],
            "kofam_file": str(kofam_path.relative_to(ROOT)) if kofam_path.exists() else "",
            "has_cbb_rbcLS_prk": yesno(CBBCORE <= present),
            "has_core_soxXABYZ": yesno(SOX_CORE <= present),
        }
        for marker in MARKERS:
            out[f"has_{marker}"] = yesno(marker in present)
            out[f"n_{marker}"] = len(by_marker[marker])
        summary_rows.append(out)

    summary_fields = [
        "genome",
        "gtdb_species",
        "omdb_completeness_pct",
        "n_proteins",
        "kofam_file",
        "has_cbb_rbcLS_prk",
        "has_core_soxXABYZ",
    ]
    for marker in MARKERS:
        summary_fields.extend([f"has_{marker}", f"n_{marker}"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    hit_fields = ["genome", "gtdb_species", "marker", "ko", "query", "raw_line"]
    with HITS_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=hit_fields)
        writer.writeheader()
        writer.writerows(hit_rows)

    print(f"Wrote {OUT.relative_to(ROOT)}")
    print(f"Wrote {HITS_OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
