"""Companion to plot_vibrio_depth_clusters.py: does the same
depth-vs-phaC-cluster pattern show up in other genera, or is it a
Vibrio-specific quirk? Three genera on the x-axis (Vibrio, plus DTSX01
and REDSEA-S09-B13 -- picked for having both a wide depth range and a
manageable number of distinct phaC clusters, same as Vibrio), depth_m on
y, points colored by cluster (jittered on x within each genus's column
so overlapping same-depth genomes are still visible).

Each genus gets its own independent color palette -- a color repeating
across genus columns is not meaningful (e.g. "orange" in the Vibrio
column and "orange" in the DTSX01 column are unrelated clusters), so the
legend is split into three per-genus blocks rather than one shared key.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from phaatlas.pipeline import sequence_clustering as sc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

GENERA = ['Vibrio', 'DTSX01', 'REDSEA-S09-B13']
CAT_PALETTE = ['#1B4F9C', '#E0621A', '#1B8A3E', '#C21F6E', '#6A3D9A', '#C2A83E',
               '#1E6E7A', '#9E3B3B', '#5B7C99', '#A8763E', '#3E9E8C', '#556B4F', '#B33951']

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
rows = list(csv.DictReader(open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline=''), delimiter='\t'))

data = {g: [] for g in GENERA}
seen = {g: set() for g in GENERA}
for r in rows:
    g = r.get('gtdb_genus', '')
    if g not in GENERA or not r.get('depth_m', ''):
        continue
    genome = r['genome']
    if genome in seen[g]:
        continue
    seen[g].add(genome)
    cid = assignments.get(r.get('target_id', ''))
    if cid:
        data[g].append((float(r['depth_m']), cid[-12:]))

colors_by_genus = {}
cluster_order_by_genus = {}
for g in GENERA:
    clusters = sorted(set(c for _, c in data[g]), key=lambda c: -sum(1 for d, cc in data[g] if cc == c))
    cluster_order_by_genus[g] = clusters
    colors_by_genus[g] = {c: CAT_PALETTE[i % len(CAT_PALETTE)] for i, c in enumerate(clusters)}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig = plt.figure(figsize=(15, 8.5), dpi=300)
ax = fig.add_axes([0.08, 0.13, 0.62, 0.68])

x_of_genus = {g: i for i, g in enumerate(GENERA)}
JITTER = 0.34

# Real depth_m values are heavily rounded (many genomes literally share one
# reported depth), so many points land exactly on top of each other -- pure
# random jitter still leaves them overlapping in clumps. Instead, group by
# exact depth and space that group's points evenly across the jitter width,
# sorted by cluster so same-color points sit next to each other (makes "how
# many of this color at this depth" readable at a glance instead of a blob).
for g in GENERA:
    x0 = x_of_genus[g]
    by_depth = {}
    for depth, cid in data[g]:
        by_depth.setdefault(depth, []).append(cid)
    for depth, cids in by_depth.items():
        cids_sorted = sorted(cids)
        n = len(cids_sorted)
        for k, cid in enumerate(cids_sorted):
            offset = 0.0 if n == 1 else -JITTER + 2 * JITTER * k / (n - 1)
            ax.scatter(x0 + offset, depth, s=42, color=colors_by_genus[g][cid], edgecolor='white',
                       linewidth=0.5, alpha=0.9, zorder=3)

ax.set_xlim(-0.6, len(GENERA) - 0.4)
ax.set_xticks(range(len(GENERA)))
ax.set_xticklabels([f'{g}\n(n={len(data[g])}, {len(cluster_order_by_genus[g])} clusters)' for g in GENERA], fontsize=10.5, fontweight='bold')
ax.set_ylabel('Depth (m)', fontsize=12)
ax.invert_yaxis()  # shallow at top, deep at bottom -- reads like a real water column
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
for i in range(len(GENERA) - 1):
    ax.axvline(i + 0.5, color='#D8DED7', linewidth=0.8, zorder=1)

fig.suptitle('Does the Depth / phaC-Cluster Pattern Hold Beyond Vibrio?', fontsize=16.5, fontweight='bold', y=0.975)
fig.text(0.5, 0.89, 'Each point is one genome; color = its phaC cluster (own palette per genus — colors are not comparable across columns)',
          ha='center', fontsize=10, color='#5B6E70')

# ---- per-genus legend blocks, stacked to the right of the plot ----
legend_x = 0.755
block_top = 0.83
for g in GENERA:
    fig.text(legend_x, block_top, g, fontsize=11.5, fontweight='bold', style='italic', color='#20302C', va='top')
    clusters = cluster_order_by_genus[g]
    n_show = min(len(clusters), 8)
    for i, cid in enumerate(clusters[:n_show]):
        y = block_top - 0.032 - i * 0.026
        n_genomes = sum(1 for d, cc in data[g] if cc == cid)
        fig.patches.append(plt.Rectangle((legend_x, y - 0.009), 0.014, 0.017, transform=fig.transFigure,
                                          facecolor=colors_by_genus[g][cid], edgecolor='#20302C', linewidth=0.5, zorder=10))
        fig.text(legend_x + 0.022, y, f'…{cid}  (n={n_genomes})', fontsize=8.3, family='monospace', va='center', color='#3A4442')
    extra = len(clusters) - n_show
    next_top_offset = 0.032 + n_show * 0.026 + (0.022 if extra > 0 else 0.012)
    if extra > 0:
        fig.text(legend_x, block_top - 0.032 - n_show * 0.026, f'+ {extra} more', fontsize=8, color='#8B9B98', style='italic', va='center')
    block_top -= next_top_offset + 0.028

footnote = (
    "Genera chosen for having both a wide depth range and a manageable number of distinct phaC clusters, like Vibrio. If clusters\n"
    "cluster tightly by depth band within a column (as in the standalone Vibrio figure), that genus shows the same depth-structuring pattern."
)
fig.text(0.5, 0.02, footnote, ha='center', va='bottom', fontsize=8.3, color='#5B6E70')

out_path = OUT / 'three_genera_depth_clusters.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'three_genera_depth_clusters.pdf', facecolor='white')
print('saved pdf too')
