"""Plot phaC genome depth against WOA23 annual temperature.

Main panel: all phaC genomes with WOA23 annual temperature and non-negative
depth. Zoom panel: subzero WOA matches, colored/labeled by common genera.
"""
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent.parent.parent
IN = ROOT / "data" / "temstapro_inputs" / "phaC_genomes_woa23_annual_temperature.tsv"
OUT = ROOT / "figures"

TOP_N = 6
ZOOM_MAX_DEPTH_M = 250
PALETTE = [
    "#1E6E7A", "#C9622D", "#8B5FBF", "#3E8914", "#C2A83E",
    "#B33951", "#4A7FB5", "#D98E04", "#6B4226", "#6FA88A",
]


def load_rows():
    rows = []
    skipped_negative_depth = 0
    with IN.open(newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r.get("woa23_temp_annual_degC") or not r.get("depth_m"):
                continue
            depth = float(r["depth_m"])
            if depth < 0:
                skipped_negative_depth += 1
                continue
            temp = float(r["woa23_temp_annual_degC"])
            genus = r.get("gtdb_genus") or "Unknown"
            rows.append({**r, "temp": temp, "depth": depth, "genus": genus})
    return rows, skipped_negative_depth


def median(vals):
    vals = sorted(vals)
    n = len(vals)
    if n == 0:
        return math.nan
    if n % 2:
        return vals[n // 2]
    return 0.5 * (vals[n // 2 - 1] + vals[n // 2])


def write_cold_summary(cold_zoom, top_genera):
    out = OUT / "phaC_woa_subzero_top_genera.tsv"
    with out.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "genus", "n_genomes", "median_woa23_temp_C", "median_depth_m",
            "min_woa23_temp_C", "max_depth_m",
        ])
        for genus in top_genera:
            g = [r for r in cold_zoom if r["genus"] == genus]
            w.writerow([
                genus, len(g), f"{median([r['temp'] for r in g]):.3f}",
                f"{median([r['depth'] for r in g]):.1f}",
                f"{min(r['temp'] for r in g):.3f}",
                f"{max(r['depth'] for r in g):.1f}",
            ])
    print(f"wrote {out}")


def main():
    rows, skipped_negative_depth = load_rows()
    cold = [r for r in rows if r["temp"] < 0]
    cold_zoom = [r for r in cold if r["depth"] <= ZOOM_MAX_DEPTH_M]
    top_genera = [g for g, _ in Counter(r["genus"] for r in cold_zoom).most_common(TOP_N)]
    color_for = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(top_genera)}
    write_cold_summary(cold_zoom, top_genera)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(13.5, 7.2), dpi=300)
    ax = fig.add_axes([0.075, 0.14, 0.55, 0.74])
    ax_zoom = fig.add_axes([0.70, 0.18, 0.26, 0.62])

    # Main panel: all WOA-matched phaC genomes.
    warm = [r for r in rows if r["temp"] >= 0]
    ax.scatter(
        [r["temp"] for r in warm], [r["depth"] for r in warm],
        s=8, c="#8FA1A0", alpha=0.20, edgecolors="none", rasterized=True,
        label="WOA >= 0 deg C",
    )
    ax.scatter(
        [r["temp"] for r in cold], [r["depth"] for r in cold],
        s=13, c="#224B5A", alpha=0.55, edgecolors="none", rasterized=True,
        label="WOA < 0 deg C",
    )
    ax.axvline(0, color="#293C40", linewidth=1.0, linestyle=(0, (4, 3)))
    ax.set_xlabel("WOA23 annual temperature at genome depth (deg C)")
    ax.set_ylabel("Sample depth (m)")
    ax.set_ylim(5650, -80)
    ax.set_xlim(-2.2, 31)
    ax.grid(True, color="#D9E0DE", linewidth=0.6, alpha=0.8)
    ax.set_title("phaC genomes across ocean temperature-depth space", fontsize=13, fontweight="bold", pad=10)

    # Zoom panel: subzero region, color top genera and gray the rest.
    other = [r for r in cold_zoom if r["genus"] not in color_for]
    ax_zoom.scatter(
        [r["temp"] for r in other], [r["depth"] for r in other],
        s=12, c="#B9C2C0", alpha=0.42, edgecolors="none", rasterized=True,
        label="Other genera",
    )
    for genus in top_genera:
        g = [r for r in cold_zoom if r["genus"] == genus]
        ax_zoom.scatter(
            [r["temp"] for r in g], [r["depth"] for r in g],
            s=24, c=color_for[genus], alpha=0.76, edgecolors="white",
            linewidths=0.25, rasterized=True,
        )

    ax_zoom.axvline(0, color="#293C40", linewidth=1.0, linestyle=(0, (4, 3)))
    ax_zoom.set_xlim(-1.55, 0.06)
    ax_zoom.set_ylim(260, -15)
    ax_zoom.grid(True, color="#D9E0DE", linewidth=0.6, alpha=0.8)
    ax_zoom.set_xlabel("WOA23 temperature (deg C)")
    ax_zoom.set_ylabel("Depth (m)")
    ax_zoom.set_title("Zoom: WOA < 0 deg C, upper 250 m", fontsize=12, fontweight="bold", pad=8)

    # Label top genera at robust central positions.
    offsets = {
        "HTCC2207": (-0.20, 34),
        "Planktomarina": (0.08, 56),
        "Marimicrobium": (-0.32, -18),
        "ASP10-02a": (0.06, -42),
        "Yoonia": (-0.27, 74),
        "Ascidiaceihabitans": (0.05, 86),
    }
    grouped = defaultdict(list)
    for r in cold_zoom:
        if r["genus"] in color_for:
            grouped[r["genus"]].append(r)
    for genus, g in grouped.items():
        x = median([r["temp"] for r in g])
        y = median([r["depth"] for r in g])
        dx, dy = offsets.get(genus, (0.05, 50))
        ax_zoom.annotate(
            f"{genus} ({len(g)})",
            xy=(x, y), xytext=(x + dx, y + dy),
            fontsize=8.4, fontstyle="italic", fontweight="bold",
            color=color_for[genus],
            arrowprops=dict(arrowstyle="-", color=color_for[genus], linewidth=0.8, alpha=0.9),
            path_effects=[pe.withStroke(linewidth=2.8, foreground="white")],
        )

    # Visual cue connecting main subzero band to zoom.
    ax.add_patch(plt.Rectangle(
        (-1.55, -15), 1.61, 275, fill=False, linewidth=1.2,
        linestyle=(0, (3, 2)), edgecolor="#293C40", alpha=0.85,
    ))

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#8FA1A0", alpha=0.55, markersize=6, label=f"All WOA-matched phaC genomes (n={len(rows):,})"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#224B5A", alpha=0.8, markersize=6, label=f"Subzero subset (n={len(cold):,})"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)

    fig.suptitle("PhaC Genomes Matched to WOA23 Annual Temperature", fontsize=16, fontweight="bold", y=0.975)
    fig.text(
        0.075, 0.045,
        "WOA23 annual climatological mean temperature; nearest 1-degree cell, vertically interpolated to sample depth.\n"
        f"Inset labels the most common genera among subzero upper-250 m matches (n={len(cold_zoom):,}); "
        f"negative-depth rows excluded (n={skipped_negative_depth}).",
        fontsize=8.6, color="#546765",
    )

    png = OUT / "phaC_woa_temperature_depth.png"
    pdf = OUT / "phaC_woa_temperature_depth.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    print(f"rows={len(rows)} cold={len(cold)} skipped_negative_depth={skipped_negative_depth}")


if __name__ == "__main__":
    main()
