#!/usr/bin/env python3
"""Query OMDB per-genome KEGG files for methane-oxidation marker genes,
for the 20 phaC-positive Methylocystis/Methylosinus genomes (see
build_methylocystis_phac_breakdown.py for how that set was derived).
Genus-level GTDB taxonomy already implies these are methanotrophs, but
this checks for the actual marker genes rather than resting on the name,
mirroring query_thioglobus_sup05_kegg_markers.py's approach for the
Thioglobus/SUP05 sulfur/carbon markers.

Marker set: pmoABC (particulate methane monooxygenase -- the core marker
nearly all aerobic methanotrophs carry), mmoXYZBCD (soluble methane
monooxygenase -- only some strains, famously including Methylosinus
trichosporium OB3b under copper limitation, one of the genomes in this
set), and mxaF/xoxF (methanol dehydrogenase, the next step after
methane oxidation to methanol).
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BREAKDOWN = ROOT / "figures" / "methylocystis_phac_taxonomic_breakdown.tsv"
CACHE = ROOT / "data" / "omdb_kegg" / "methylocystis_methanotrophs"
GENOME_PAGE_CACHE = ROOT / "data" / "omdb_genome_pages" / "methylocystis_methanotrophs"
PRESENCE_OUT = ROOT / "figures" / "methylocystis_methanotrophy_marker_presence.tsv"
HITS_OUT = ROOT / "figures" / "methylocystis_methanotrophy_marker_hits.tsv"

MARKERS = {
    "pmoA": {"K10944"},
    "pmoB": {"K10945"},
    "pmoC": {"K10946"},
    "mmoX": {"K16157"},
    "mmoY": {"K16158"},
    "mmoZ": {"K16159"},
    "mmoB": {"K16161"},
    "mmoC": {"K16162"},
    "mmoD": {"K16163"},
    "mxaF": {"K14028"},
    "xoxF": {"K23995"},
}
PMO_CORE = {"pmoA", "pmoB", "pmoC"}
MMO_CORE = {"mmoX", "mmoY", "mmoZ"}
GENOME_JSON_RE = re.compile(r"window\._MB_GENOME_DATA = JSON\.parse\(`(.+?)`\);", re.S)


def api_url(genome: str) -> str:
    return (
        "https://motus-api.microbiomics.io/v1/genomes/"
        f"{urllib.parse.quote(genome)}"
        "/download?file_type=kegg"
    )


def genome_page_url(genome: str) -> str:
    return f"https://omdb.microbiomics.io/repository/ocean/genome/{urllib.parse.quote(genome)}"


def fetch_url(url: str, timeout: int = 90) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "pha-ocean-atlas-marker-query/1.0"})
    context = ssl._create_unverified_context()
    with urllib.request.urlopen(req, timeout=timeout, context=context) as resp:
        return resp.read()


def download_kegg(genome: str, out: Path) -> str:
    if out.exists() and out.stat().st_size > 0:
        return "cached"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    try:
        tmp.write_bytes(fetch_url(api_url(genome)))
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


def fetch_genome_metadata(genome: str) -> tuple[dict[str, object], str]:
    out = GENOME_PAGE_CACHE / f"{genome}.html"
    status = "cached"
    if not out.exists() or out.stat().st_size == 0:
        out.parent.mkdir(parents=True, exist_ok=True)
        status = "downloaded"
        try:
            out.write_bytes(fetch_url(genome_page_url(genome)))
        except urllib.error.HTTPError as exc:
            return {}, f"http_{exc.code}"
        except urllib.error.URLError as exc:
            return {}, f"url_error:{exc.reason}"
        time.sleep(0.15)
    text = out.read_text(errors="replace")
    match = GENOME_JSON_RE.search(text)
    if not match:
        return {}, "json_not_found"
    return json.loads(match.group(1)), status


def read_kegg_hits(path: Path) -> list[dict[str, str]]:
    hits = []
    with gzip.open(path, "rt", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ko = row.get("KO", "")
            if ko == "NA":
                continue
            for marker, kos in MARKERS.items():
                if ko in kos:
                    hits.append({
                        "query": row.get("QUERY", ""), "marker": marker, "ko": ko,
                        "description": row.get("DESCRIPTION", ""),
                    })
    return hits


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def fmt_float(value: object) -> str:
    if value in (None, ""):
        return ""
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value)


def main() -> int:
    with BREAKDOWN.open() as handle:
        genomes = list(csv.DictReader(handle, delimiter="\t"))

    presence_rows = []
    hit_rows = []
    for row in genomes:
        genome = row["genome"]
        kegg_path = CACHE / f"{genome}.kegg.apr22.tsv.gz"
        status = download_kegg(genome, kegg_path)
        genome_meta, genome_meta_status = fetch_genome_metadata(genome)
        hits = read_kegg_hits(kegg_path) if kegg_path.exists() and status in {"cached", "downloaded"} else []
        marker_counts = Counter(hit["marker"] for hit in hits)
        present = {marker for marker, count in marker_counts.items() if count > 0}

        for hit in hits:
            hit_rows.append({
                "genome": genome, "gtdb_genus": row["gtdb_genus"], "gtdb_species": row["gtdb_species"],
                **hit,
            })

        presence_rows.append({
            "genome": genome,
            "gtdb_genus": row["gtdb_genus"],
            "gtdb_species": row["gtdb_species"],
            "n_qc_targets": row["n_qc_targets"],
            "any_triad_complete_alignment": row["any_triad_complete_alignment"],
            "omdb_completeness_pct": fmt_float(genome_meta.get("completeness")),
            "omdb_contamination_pct": fmt_float(genome_meta.get("contamination")),
            "has_pmoA": yesno("pmoA" in present),
            "has_pmoB": yesno("pmoB" in present),
            "has_pmoC": yesno("pmoC" in present),
            "has_pmoABC_core": yesno(PMO_CORE <= present),
            "has_mmoX": yesno("mmoX" in present),
            "has_mmoY": yesno("mmoY" in present),
            "has_mmoZ": yesno("mmoZ" in present),
            "has_mmoBCD": yesno({"mmoB", "mmoC", "mmoD"} & present != set()),
            "has_mmoXYZ_core": yesno(MMO_CORE <= present),
            "has_mxaF": yesno("mxaF" in present),
            "has_xoxF": yesno("xoxF" in present),
            "markers_present": ",".join(sorted(present)),
            "kegg_download_status": status,
            "genome_metadata_status": genome_meta_status,
        })
        print(f"{genome}: {sorted(present) or 'no markers found'} (kegg={status}, meta={genome_meta_status})")

    presence_fields = list(presence_rows[0].keys()) if presence_rows else []
    PRESENCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    with PRESENCE_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=presence_fields)
        writer.writeheader()
        writer.writerows(presence_rows)

    hit_fields = ["genome", "gtdb_genus", "gtdb_species", "query", "marker", "ko", "description"]
    with HITS_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=hit_fields)
        writer.writeheader()
        writer.writerows(hit_rows)

    n_pmo = sum(1 for r in presence_rows if r["has_pmoABC_core"] == "yes")
    n_mmo = sum(1 for r in presence_rows if r["has_mmoXYZ_core"] == "yes")
    n_any = sum(1 for r in presence_rows if r["markers_present"])
    print()
    print(f"{n_any}/{len(presence_rows)} genomes have >=1 methanotrophy marker")
    print(f"{n_pmo}/{len(presence_rows)} genomes have the full pmoABC core (particulate MMO)")
    print(f"{n_mmo}/{len(presence_rows)} genomes have the full mmoXYZ core (soluble MMO)")
    print(f"Wrote {PRESENCE_OUT}")
    print(f"Wrote {HITS_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
