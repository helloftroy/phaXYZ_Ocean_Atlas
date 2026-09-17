"""Plot % of all OMDB genomes carrying phaC by geographic grid cell.

Denominator: all OMDB genomes with sample latitude/longitude.
Numerator: distinct genomes resolved from phaC NR100 clusters. As a
fallback, a phaC *_with_metadata* table can be used, but that may be
capped depending on how metadata enrichment was run.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_ALL = ROOT / "data" / "all_genomes" / "omdb_all_genomes_with_locations.tsv"
DEFAULT_PHAC_GENOMES = ROOT / "data" / "all_genomes" / "phaC_all_genomes_from_nr100_clusters.tsv"
DEFAULT_PHAC = ROOT / "PHA_bioprospecting" / "omdb_search" / "results" / "phaC_unique_targets_with_metadata_depth.tsv"
DEFAULT_OUT_TSV = ROOT / "data" / "all_genomes" / "phaC_prevalence_all_genomes_by_region.tsv"
DEFAULT_OUT_PREFIX = ROOT / "figures" / "phaC_pct_all_genomes_by_region"


def parse_float(raw: str):
    if raw is None or raw == "":
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def cell_for(lat: float, lon: float, grid_deg: int) -> tuple[int, int]:
    lat_bin = math.floor(lat / grid_deg) * grid_deg
    lon_bin = math.floor(lon / grid_deg) * grid_deg
    return lat_bin, lon_bin


def load_genomes_from_tsv(path: Path) -> set[str]:
    genomes = set()
    with path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genome = row.get("genome", "")
            if genome:
                genomes.add(genome)
    return genomes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-genomes", type=Path, default=DEFAULT_ALL)
    ap.add_argument("--phac-genomes", type=Path, default=DEFAULT_PHAC_GENOMES)
    ap.add_argument("--phac-metadata", type=Path, default=DEFAULT_PHAC)
    ap.add_argument("--out-tsv", type=Path, default=DEFAULT_OUT_TSV)
    ap.add_argument("--out-prefix", type=Path, default=DEFAULT_OUT_PREFIX)
    ap.add_argument("--grid-deg", type=int, default=10)
    ap.add_argument("--min-genomes-per-cell", type=int, default=20)
    args = ap.parse_args()

    if not args.all_genomes.exists():
        raise SystemExit(
            f"All-genome location table not found: {args.all_genomes}\n"
            "Build it first with figures/scripts/export_omdb_all_genome_locations.py"
        )
    if args.phac_genomes.exists():
        phac_source = args.phac_genomes
    elif args.phac_metadata.exists():
        phac_source = args.phac_metadata
        print(
            "WARNING: using phaC metadata table as numerator fallback. "
            "For an uncapped numerator, first run "
            "figures/scripts/export_phac_genomes_from_omdb_clusters.py"
        )
    else:
        raise SystemExit(
            f"Neither phaC genome table nor metadata table found:\n"
            f"  {args.phac_genomes}\n"
            f"  {args.phac_metadata}"
        )

    phac_genomes = load_genomes_from_tsv(phac_source)
    cell_total = defaultdict(int)
    cell_phac = defaultdict(int)
    total_genomes_with_location = 0
    phac_genomes_with_location = 0
    seen_genomes = set()

    with args.all_genomes.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genome = row.get("genome", "")
            if not genome or genome in seen_genomes:
                continue
            seen_genomes.add(genome)
            lat = parse_float(row.get("latitude_degN", ""))
            lon = parse_float(row.get("longitude_degE", ""))
            if lat is None or lon is None:
                continue
            cell = cell_for(lat, lon, args.grid_deg)
            cell_total[cell] += 1
            total_genomes_with_location += 1
            if genome in phac_genomes:
                cell_phac[cell] += 1
                phac_genomes_with_location += 1

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "lat_bin", "lon_bin", "n_all_genomes", "n_phac_genomes",
            "pct_all_genomes_with_phac", "passes_min_genomes_filter",
            "phac_numerator_source",
        ])
        for cell in sorted(cell_total):
            total = cell_total[cell]
            phac = cell_phac[cell]
            writer.writerow([
                cell[0], cell[1], total, phac, f"{100 * phac / total:.6f}",
                "yes" if total >= args.min_genomes_per_cell else "no",
                phac_source,
            ])

    lat_bins = list(range(-90, 90, args.grid_deg))
    lon_bins = list(range(-180, 180, args.grid_deg))
    lat_idx = {b: i for i, b in enumerate(lat_bins)}
    lon_idx = {b: i for i, b in enumerate(lon_bins)}
    pct_grid = np.full((len(lat_bins), len(lon_bins)), np.nan)
    n_grid = np.zeros((len(lat_bins), len(lon_bins)))
    for (lat_bin, lon_bin), total in cell_total.items():
        if lat_bin not in lat_idx or lon_bin not in lon_idx:
            continue
        i, j = lat_idx[lat_bin], lon_idx[lon_bin]
        n_grid[i, j] = total
        if total >= args.min_genomes_per_cell:
            pct_grid[i, j] = 100 * cell_phac[(lat_bin, lon_bin)] / total

    observed = pct_grid[~np.isnan(pct_grid)]
    if observed.size == 0:
        raise SystemExit("No grid cells passed the minimum-genome filter.")

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "phac_all_genomes_pct",
        ["#F5F1E8", "#CAD8D1", "#77A6A7", "#1E6E7A", "#153F4A"],
    )
    lon_edges = lon_bins + [180]
    lat_edges = lat_bins + [90]
    masked = np.ma.masked_invalid(pct_grid)
    insufficient = (n_grid > 0) & np.isnan(pct_grid)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(figsize=(15, 7.5), dpi=300)
    ax.set_facecolor("#F3F5F2")
    mesh = ax.pcolormesh(
        lon_edges, lat_edges, masked,
        cmap=cmap, vmin=0, vmax=float(np.nanpercentile(observed, 98)),
        edgecolors="white", linewidth=0.35,
    )
    for i, lat_bin in enumerate(lat_bins):
        for j, lon_bin in enumerate(lon_bins):
            if insufficient[i, j]:
                ax.add_patch(plt.Rectangle(
                    (lon_bin, lat_bin), args.grid_deg, args.grid_deg,
                    facecolor="#D8DED7", hatch="////", edgecolor="#B9C3BB",
                    linewidth=0.2,
                ))

    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 85)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="#BCC8C3", linewidth=0.45, linestyle=(0, (1, 4)))
    ax.set_title(
        "phaC Prevalence Among All OMDB Genomes",
        fontsize=17, fontweight="bold", pad=14,
    )
    fig.text(
        0.5, 0.90,
        f"Distinct phaC-positive genomes divided by all OMDB genomes with sample lat/lon, "
        f"{args.grid_deg}x{args.grid_deg} degree cells; hatched cells have <{args.min_genomes_per_cell} genomes.",
        ha="center", fontsize=10.5, color="#566A67",
    )
    cbar = fig.colorbar(mesh, ax=ax, orientation="vertical", shrink=0.78, pad=0.018)
    cbar.set_label("% of all genomes with phaC")
    fig.text(
        0.5, 0.035,
        f"All genomes with location: {total_genomes_with_location:,}; "
        f"phaC-positive among those: {phac_genomes_with_location:,}; "
        f"global prevalence: {100 * phac_genomes_with_location / total_genomes_with_location:.2f}%. "
        f"Color scale capped at the 98th percentile of shown cells ({np.nanpercentile(observed, 98):.2f}%).\n"
        f"Numerator source: {phac_source.name}",
        ha="center", fontsize=8.8, color="#566A67",
    )

    png = args.out_prefix.with_suffix(".png")
    pdf = args.out_prefix.with_suffix(".pdf")
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    print(f"wrote {args.out_tsv}")
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    print(f"phaC numerator source: {phac_source}")
    print(f"all genomes with location: {total_genomes_with_location:,}")
    print(f"phaC genomes with location: {phac_genomes_with_location:,}")
    print(f"global %: {100 * phac_genomes_with_location / total_genomes_with_location:.3f}")
    print(f"cells shown: {int(np.sum(~np.isnan(pct_grid))):,}")


if __name__ == "__main__":
    main()
