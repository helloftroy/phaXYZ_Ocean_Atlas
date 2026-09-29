#!/usr/bin/env python3
"""Query OMDB per-genome KEGG files for Thioglobus/SUP05 marker genes."""

from __future__ import annotations

import csv
import gzip
import json
import re
import ssl
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BREAKDOWN = ROOT / "figures" / "thioglobus_sup05_phac_taxonomic_breakdown.tsv"
CACHE = ROOT / "data" / "omdb_kegg" / "thioglobus_sup05"
GENOME_PAGE_CACHE = ROOT / "data" / "omdb_genome_pages" / "thioglobus_sup05"
PRESENCE_OUT = ROOT / "figures" / "thioglobus_sup05_kegg_marker_presence.tsv"
HITS_OUT = ROOT / "figures" / "thioglobus_sup05_kegg_marker_hits.tsv"
SUMMARY_OUT = ROOT / "figures" / "thioglobus_sup05_kegg_marker_summary_by_subclade.tsv"
SPECIES_SUMMARY_OUT = ROOT / "figures" / "thioglobus_sup05_kegg_marker_summary_by_species.tsv"


MARKERS = {
    "rbcL": {"K01601"},
    "rbcS": {"K01602"},
    "prk": {"K00855"},
    "sqr": {"K17218"},
    "fccB": {"K17229"},
    "fccA": {"K17230"},
    "soxA": {"K17222"},
    "soxX": {"K17223"},
    "soxB": {"K17224"},
    "soxC": {"K17225"},
    "soxD": {"K22622"},
    "soxY": {"K17226"},
    "soxZ": {"K17227"},
}

SOX_CORE = {"soxX", "soxA", "soxB", "soxY", "soxZ"}
SOX_WITH_CD = SOX_CORE | {"soxC", "soxD"}
GENOME_JSON_RE = re.compile(r"window\._MB_GENOME_DATA = JSON\.parse\(`(.+?)`\);", re.S)


def subclade(row: dict[str, str]) -> str:
    genus = row["gtdb_genus"]
    if genus == "Thioglobus_A":
        return "Thioglobus_A chemoautotrophic"
    if genus in {"Thioglobus", "Thioglobus_B"}:
        return "Thioglobus_B / Thioglobus sensu stricto mixotrophic"
    if genus == "Pseudothioglobus":
        return "Pseudothioglobus mixotrophic-to-heterotrophic"
    return genus or "unknown"


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
                    hits.append(
                        {
                            "query": row.get("QUERY", ""),
                            "marker": marker,
                            "ko": ko,
                            "description": row.get("DESCRIPTION", ""),
                            "module": row.get("MODULE", ""),
                            "pathway": row.get("PATHWAY", ""),
                        }
                    )
    return hits


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def fmt_float(value: object) -> str:
    if value in (None, ""):
        return ""
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value)


def mean_float(rows: list[dict[str, object]], key: str) -> str:
    vals = []
    for row in rows:
        try:
            vals.append(float(row[key]))
        except (TypeError, ValueError):
            pass
    return f"{sum(vals) / len(vals):.2f}" if vals else ""


def main() -> int:
    with BREAKDOWN.open() as handle:
        genomes = list(csv.DictReader(handle, delimiter="\t"))

    presence_rows: list[dict[str, object]] = []
    hit_rows: list[dict[str, object]] = []

    for row in genomes:
        genome = row["genome"]
        kegg_path = CACHE / f"{genome}.kegg.apr22.tsv.gz"
        status = download_kegg(genome, kegg_path)
        genome_meta, genome_meta_status = fetch_genome_metadata(genome)
        hits = read_kegg_hits(kegg_path) if kegg_path.exists() and status in {"cached", "downloaded"} else []
        marker_counts = Counter(hit["marker"] for hit in hits)
        present = {marker for marker, count in marker_counts.items() if count > 0}
        sox_present = sorted(m for m in present if m.startswith("sox"))

        for hit in hits:
            hit_rows.append(
                {
                    "genome": genome,
                    "gtdb_genus": row["gtdb_genus"],
                    "gtdb_species": row["gtdb_species"],
                    "ecological_subclade": subclade(row),
                    "query": hit["query"],
                    "marker": hit["marker"],
                    "ko": hit["ko"],
                    "description": hit["description"],
                    "module": hit["module"],
                    "pathway": hit["pathway"],
                }
            )

        cbb_markers = [m for m in ("rbcL", "rbcS", "prk") if m in present]
        fcc_markers = [m for m in ("fccA", "fccB") if m in present]
        presence_rows.append(
            {
                "genome": genome,
                "gtdb_genus": row["gtdb_genus"],
                "gtdb_species": row["gtdb_species"],
                "ecological_subclade": subclade(row),
                "omdb_completeness_pct": fmt_float(genome_meta.get("completeness")),
                "omdb_contamination_pct": fmt_float(genome_meta.get("contamination")),
                "checkm_completeness_pct": fmt_float(genome_meta.get("checkm_completeness")),
                "anvio_completion_pct": fmt_float(genome_meta.get("anvio_completion")),
                "n_phac_targets": row["n_phac_targets"],
                "n_triad_complete": row["n_triad_complete"],
                "has_rbcL": yesno("rbcL" in present),
                "has_rbcS": yesno("rbcS" in present),
                "has_prk": yesno("prk" in present),
                "cbb_markers": ",".join(cbb_markers),
                "has_cbb_rbcLS_prk": yesno(set(cbb_markers) == {"rbcL", "rbcS", "prk"}),
                "has_sqr": yesno("sqr" in present),
                "has_fccA": yesno("fccA" in present),
                "has_fccB": yesno("fccB" in present),
                "has_fccAB": yesno(set(fcc_markers) == {"fccA", "fccB"}),
                "has_soxX": yesno("soxX" in present),
                "has_soxA": yesno("soxA" in present),
                "has_soxB": yesno("soxB" in present),
                "has_soxY": yesno("soxY" in present),
                "has_soxZ": yesno("soxZ" in present),
                "has_soxC": yesno("soxC" in present),
                "has_soxD": yesno("soxD" in present),
                "sox_markers": ",".join(sox_present),
                "n_sox_markers": len(sox_present),
                "has_core_soxXABYZ": yesno(SOX_CORE <= present),
                "has_sox_module_with_cd": yesno(SOX_WITH_CD <= present),
                "kegg_download_status": status,
                "kegg_file": str(kegg_path.relative_to(ROOT)) if kegg_path.exists() else "",
                "genome_metadata_status": genome_meta_status,
            }
        )

    presence_fields = [
        "genome",
        "gtdb_genus",
        "gtdb_species",
        "ecological_subclade",
        "omdb_completeness_pct",
        "omdb_contamination_pct",
        "checkm_completeness_pct",
        "anvio_completion_pct",
        "n_phac_targets",
        "n_triad_complete",
        "has_rbcL",
        "has_rbcS",
        "has_prk",
        "cbb_markers",
        "has_cbb_rbcLS_prk",
        "has_sqr",
        "has_fccA",
        "has_fccB",
        "has_fccAB",
        "has_soxX",
        "has_soxA",
        "has_soxB",
        "has_soxY",
        "has_soxZ",
        "has_soxC",
        "has_soxD",
        "sox_markers",
        "n_sox_markers",
        "has_core_soxXABYZ",
        "has_sox_module_with_cd",
        "kegg_download_status",
        "kegg_file",
        "genome_metadata_status",
    ]
    hit_fields = [
        "genome",
        "gtdb_genus",
        "gtdb_species",
        "ecological_subclade",
        "query",
        "marker",
        "ko",
        "description",
        "module",
        "pathway",
    ]
    write_tsv(PRESENCE_OUT, presence_rows, presence_fields)
    write_tsv(HITS_OUT, hit_rows, hit_fields)

    by_subclade: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in presence_rows:
        by_subclade[str(row["ecological_subclade"])].append(row)

    summary_rows: list[dict[str, object]] = []
    summary_markers = [
        "has_cbb_rbcLS_prk",
        "has_rbcL",
        "has_rbcS",
        "has_prk",
        "has_sqr",
        "has_fccA",
        "has_fccB",
        "has_fccAB",
        "has_soxX",
        "has_soxA",
        "has_soxB",
        "has_soxY",
        "has_soxZ",
        "has_soxC",
        "has_soxD",
        "has_core_soxXABYZ",
        "has_sox_module_with_cd",
    ]
    for name, rows in sorted(by_subclade.items()):
        out = {"ecological_subclade": name, "n_genomes": len(rows)}
        out["mean_omdb_completeness_pct"] = mean_float(rows, "omdb_completeness_pct")
        out["min_omdb_completeness_pct"] = min(
            (float(r["omdb_completeness_pct"]) for r in rows if r["omdb_completeness_pct"]),
            default="",
        )
        if out["min_omdb_completeness_pct"] != "":
            out["min_omdb_completeness_pct"] = f"{out['min_omdb_completeness_pct']:.2f}"
        out["n_triad_complete_phac_targets"] = sum(int(r["n_triad_complete"]) for r in rows)
        for marker in summary_markers:
            out[f"n_{marker}"] = sum(1 for r in rows if r[marker] == "yes")
        summary_rows.append(out)

    summary_fields = [
        "ecological_subclade",
        "n_genomes",
        "mean_omdb_completeness_pct",
        "min_omdb_completeness_pct",
        "n_triad_complete_phac_targets",
    ] + [
        f"n_{m}" for m in summary_markers
    ]
    write_tsv(SUMMARY_OUT, summary_rows, summary_fields)

    by_species: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in presence_rows:
        by_species[
            (
                str(row["gtdb_genus"]),
                str(row["gtdb_species"]),
                str(row["ecological_subclade"]),
            )
        ].append(row)

    species_rows: list[dict[str, object]] = []
    for (genus, species, clade), rows in sorted(by_species.items()):
        out = {
            "gtdb_genus": genus,
            "gtdb_species": species,
            "ecological_subclade": clade,
            "n_genomes": len(rows),
            "mean_omdb_completeness_pct": mean_float(rows, "omdb_completeness_pct"),
            "min_omdb_completeness_pct": min(
                (float(r["omdb_completeness_pct"]) for r in rows if r["omdb_completeness_pct"]),
                default="",
            ),
            "n_triad_complete_phac_targets": sum(int(r["n_triad_complete"]) for r in rows),
        }
        if out["min_omdb_completeness_pct"] != "":
            out["min_omdb_completeness_pct"] = f"{out['min_omdb_completeness_pct']:.2f}"
        for marker in summary_markers:
            out[f"n_{marker}"] = sum(1 for r in rows if r[marker] == "yes")
        species_rows.append(out)
    species_fields = [
        "gtdb_genus",
        "gtdb_species",
        "ecological_subclade",
        "n_genomes",
        "mean_omdb_completeness_pct",
        "min_omdb_completeness_pct",
        "n_triad_complete_phac_targets",
    ] + [
        f"n_{m}" for m in summary_markers
    ]
    write_tsv(SPECIES_SUMMARY_OUT, species_rows, species_fields)

    print(f"Downloaded/cached KEGG files for {len(presence_rows)} genomes")
    print(f"Marker hits: {len(hit_rows)}")
    print(f"Wrote {PRESENCE_OUT.relative_to(ROOT)}")
    print(f"Wrote {HITS_OUT.relative_to(ROOT)}")
    print(f"Wrote {SUMMARY_OUT.relative_to(ROOT)}")
    print(f"Wrote {SPECIES_SUMMARY_OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
