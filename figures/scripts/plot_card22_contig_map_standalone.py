"""Standalone CARD22-1 Desulfoluna phaC contig map.

This redraws panel C from figures/card22_desulfoluna_phac_deep_dive.png as a
compact, panel-ready figure. It uses the already-generated
card22_desulfoluna_phac_paralogs.tsv for phaC positions/classes and the
documented one-off neighboring phaA/B/E/J positions recorded in the original
deep-dive script.

Outputs:
    figures/card22_desulfoluna_contig_map_compact.png
    figures/card22_desulfoluna_contig_map_compact.pdf
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches


ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "figures"
IN = OUT / "card22_desulfoluna_phac_paralogs.tsv"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _deep_dive_copies  # noqa: E402

OTHER_PHA_GENES_ON_SCAFFOLDS = [
    # scaffold, gene_index, family, note
    ("scaffold_1", 31, "phaB", ""),
    ("scaffold_1", 56, "phaB", ""),
    ("scaffold_1", 125, "phaB", ""),
    ("scaffold_1", 236, "phaA", ""),
    ("scaffold_1", 1050, "phaA", ""),
    ("scaffold_1", 1076, "phaA", ""),
    ("scaffold_1", 1099, "phaA", ""),
    ("scaffold_12", 81, "phaE", "upstream of phaC"),
    ("scaffold_12", 83, "phaJ", "downstream of phaC"),
]

FAMILY_COLORS = {
    "phaA": "#4A7FB5",
    "phaB": "#D98E04",
    "phaC": "#2F6F63",
    "phaE": "#7A5FA0",
    "phaJ": "#9E3B3B",
}
INK = "#20302C"
MUTED = "#667879"
GRID = "#D8E2DD"
TRACK = "#B9C7C1"
WARN = "#B33951"


GENOME = "CARD22-1_SAMN24292811_MAG_00000010"


def load_genes():
    # Only the triad-complete copies, and in the same order and colours the
    # identity heatmap and the circos use. The map previously drew all 7
    # QC-passing copies while the heatmap drew 6 and the circos 5, so the three
    # panels of one figure disagreed about how many phaC this genome has. 187467
    # (541 aa) fails the triad check by both the structural geometry and the
    # alignment columns, so it is not a functional synthase.
    styles = {c.target_id: c for c in _deep_dive_copies.copies_for(GENOME)}
    genes_by_scaffold = {}
    with IN.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["target_id"] not in styles:
                continue
            style = styles[row["target_id"]]
            scaf = row["own_scaffold"]
            idx = int(row["own_gene_index"])
            short = row["target_id"][-6:]
            triad = row["triad_complete"] == "True"
            cls = row["best_class_relative"].replace("Class", "Cl.")
            trusted = row["passes_class_trusted_cutoff"] == "True"
            label = f"phaC {short}"
            sublabel = f"{row['protein_length_aa']} aa, {cls}{'*' if trusted else ''}"
            genes_by_scaffold.setdefault(scaf, []).append({
                "idx": idx,
                "family": "phaC",
                "label": label,
                "sublabel": sublabel,
                "triad": triad,
                "ambiguous": False,
                "color": FAMILY_COLORS["phaC"],
                "chip": style.color,
            })

    for scaf, idx, fam, note in OTHER_PHA_GENES_ON_SCAFFOLDS:
        genes_by_scaffold.setdefault(scaf, []).append({
            "idx": idx,
            "family": fam,
            "label": fam,
            "sublabel": note,
            "triad": None,
            "ambiguous": False,
            "color": FAMILY_COLORS[fam],
        })

    for genes in genes_by_scaffold.values():
        genes.sort(key=lambda g: g["idx"])
    return genes_by_scaffold


def arrow(ax, x, y, direction, color, edge, scale=1.0, zorder=4):
    length = 0.78 * scale
    dx = length if direction >= 0 else -length
    start = x - dx / 2
    end = x + dx / 2
    patch = FancyArrowPatch(
        (start, y), (end, y),
        arrowstyle="Simple,tail_width=0.44,head_width=0.82,head_length=0.32",
        mutation_scale=12,
        linewidth=1.1,
        facecolor=color,
        edgecolor=edge,
        zorder=zorder,
    )
    ax.add_patch(patch)


def label_gene(ax, x, y, gene, level=0, compact=False):
    if compact and gene["family"] == "phaC":
        base_dy = 0.42
    else:
        base_dy = 0.31 if compact else 0.30
    dy = base_dy + level * (0.23 if compact else 0.27)
    weight = "bold" if gene["family"] == "phaC" else "semibold"
    fs = (7.5 if compact else 8.7) if gene["family"] == "phaC" else (7.2 if compact else 8.2)
    ax.text(
        x, y + dy, gene["label"],
        ha="center", va="bottom", fontsize=fs, fontweight=weight, color=INK,
        path_effects=[pe.withStroke(linewidth=2.5, foreground="white", alpha=0.95)],
        zorder=6,
    )
    if gene["sublabel"] and gene["family"] == "phaC":
        ax.text(
            x, y + dy - (0.19 if compact else 0.22), gene["sublabel"],
            ha="center", va="bottom", fontsize=5.8 if compact else 6.7, color=MUTED,
            path_effects=[pe.withStroke(linewidth=2.2, foreground="white", alpha=0.9)],
            zorder=6,
        )


def compressed_positions(indices, max_gap=8.0):
    """Preserve gene order while compressing long empty gaps for display."""
    ordered = sorted(indices)
    if not ordered:
        return {}
    display = {ordered[0]: 0.0}
    for prev, cur in zip(ordered, ordered[1:]):
        display[cur] = display[prev] + min(float(cur - prev), max_gap)
    return display


def draw_track(ax, genes, xlim, title, subtitle="", show_ticks=True, max_gap=8.0):
    visible = [g for g in genes if xlim[0] <= g["idx"] <= xlim[1]]
    display = compressed_positions([g["idx"] for g in visible], max_gap=max_gap)
    if display:
        disp_vals = list(display.values())
        dxlim = (min(disp_vals) - 1.2, max(disp_vals) + 1.2)
    else:
        dxlim = (0, 1)

    ax.set_xlim(*dxlim)
    ax.set_ylim(-0.58, 1.05)
    ax.axhline(0, color=TRACK, linewidth=1.1, zorder=1)
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("#A9B8B1")
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(axis="x", labelsize=7.4, colors="#465854", length=2.5, pad=2)
    if not show_ticks:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.5, alpha=0.55, zorder=0)
    ax.set_axisbelow(True)
    ax.text(0.0, 1.10, title, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=9.7, fontweight="bold", color=INK)
    if subtitle:
        ax.text(1.0, 1.10, subtitle, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.2, color=MUTED)

    if show_ticks and display:
        ordered_idxs = sorted(display)
        ax.set_xticks([display[i] for i in ordered_idxs])
        ax.set_xticklabels([str(i) for i in ordered_idxs])

    prev_x = None
    level = 0
    span = dxlim[1] - dxlim[0]
    for gene in visible:
        x = display[gene["idx"]]
        if prev_x is not None and abs(x - prev_x) < span * 0.13:
            level = (level + 1) % 3
        else:
            level = 0
        prev_x = x
        edge = INK
        arrow(ax, x, 0, 1, gene["color"], edge)
        label_gene(ax, x, 0, gene, level)


def draw_compact_track(ax, genes, title, subtitle="", max_gap=1.55, show_breaks=True):
    """Draw a short, panel-friendly contig neighborhood.

    X positions are compressed gene-index positions. Tick labels remain the
    original gene indices so the biology is still traceable in a tiny panel.
    """
    display = compressed_positions([g["idx"] for g in genes], max_gap=max_gap)
    ordered = sorted(genes, key=lambda g: g["idx"])
    xs = [display[g["idx"]] for g in ordered]
    pad = 0.65
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(-0.44, 1.08)
    ax.axis("off")

    for left, right, g0, g1 in zip(xs, xs[1:], ordered, ordered[1:]):
        mid = (left + right) / 2
        if show_breaks and g1["idx"] - g0["idx"] > 35:
            ax.plot([left + 0.36, mid - 0.16], [0, 0], color=TRACK, linewidth=1.0, zorder=1)
            ax.plot([mid + 0.16, right - 0.36], [0, 0], color=TRACK, linewidth=1.0, zorder=1)
            ax.text(mid, 0.0, "//", ha="center", va="center", fontsize=10.5,
                    fontweight="bold", color="#95A8A1", zorder=2)
        else:
            ax.plot([left + 0.36, right - 0.36], [0, 0], color=TRACK, linewidth=1.0, zorder=1)
    ax.plot([xs[0] - 0.55, xs[0] - 0.36], [0, 0], color=TRACK, linewidth=1.0, zorder=1)
    ax.plot([xs[-1] + 0.36, xs[-1] + 0.55], [0, 0], color=TRACK, linewidth=1.0, zorder=1)

    ax.text(0.0, 0.98, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=8.4, fontweight="bold", color=INK)
    if subtitle:
        ax.text(1.0, 0.98, subtitle, transform=ax.transAxes, ha="right", va="top",
                fontsize=6.3, color=MUTED)

    prev_x = None
    level = 0
    for gene in ordered:
        x = display[gene["idx"]]
        if prev_x is not None and abs(x - prev_x) < 1.75:
            level = (level + 1) % 3
        else:
            level = 0
        prev_x = x
        edge = INK
        arrow(ax, x, 0, 1, gene["color"], edge, scale=0.82)
        label_gene(ax, x, 0, gene, level, compact=True)
        ax.text(x, -0.25, str(gene["idx"]), ha="center", va="top", fontsize=5.8, color=MUTED)
        # Colour chip keying this copy to the identity heatmap's axis chips.
        # Below the arrow, not above: stacked over the label it landed on the
        # text, and the labels sit at several heights to avoid each other.
        if gene.get("chip"):
            ax.scatter([x], [-0.52], s=62, color=gene["chip"], edgecolors="white",
                       linewidths=0.9, clip_on=False, zorder=7)


def main():
    genes_by_scaffold = load_genes()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig = plt.figure(figsize=(7.0, 4.1), dpi=300)
    gs = fig.add_gridspec(
        4, 1,
        height_ratios=[1.05, 0.86, 0.86, 0.22],
        hspace=0.42,
        left=0.045,
        right=0.985,
        top=0.965,
        bottom=0.12,
    )

    scaf1 = genes_by_scaffold["scaffold_1"]
    ax = fig.add_subplot(gs[0])
    draw_compact_track(ax, scaf1, "scaffold_1", max_gap=1.42)

    small_specs = [
        ("scaffold_2", (294, 310), "single Class III-like phaC"),
        ("scaffold_9", (68, 82), "adjacent phaC pair"),
        ("scaffold_12", (76, 88), "compact phaE-phaC-phaJ cassette"),
        ("scaffold_27", (24, 40), "isolated phaC/Z-ambiguous copy"),
    ]
    small_gs = gs[1:3].subgridspec(2, 2, hspace=0.48, wspace=0.24)
    for i, (scaf, xlim, subtitle) in enumerate(small_specs):
        ax = fig.add_subplot(small_gs[i // 2, i % 2])
        visible = [g for g in genes_by_scaffold[scaf] if xlim[0] <= g["idx"] <= xlim[1]]
        draw_compact_track(ax, visible, scaf, max_gap=1.2, show_breaks=False)

    handles = [
        mpatches.Patch(facecolor=FAMILY_COLORS[f], edgecolor=INK, label=f)
        for f in ["phaA", "phaB", "phaC", "phaE", "phaJ"]
    ]
    legend_ax = fig.add_subplot(gs[3])
    legend_ax.axis("off")
    legend_ax.legend(handles=handles, loc="center", ncol=5, frameon=False,
                     fontsize=6.8, handlelength=1.0, columnspacing=0.9)

    png = OUT / "card22_desulfoluna_contig_map_compact.png"
    pdf = OUT / "card22_desulfoluna_contig_map_compact.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    print(f"saved {png}")
    print(f"saved {pdf}")


if __name__ == "__main__":
    main()
