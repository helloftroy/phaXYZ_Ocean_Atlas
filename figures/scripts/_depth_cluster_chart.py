"""Shared builder for the single-genus depth-vs-phaC-cluster chart (used by
plot_vibrio_depth_clusters.py and plot_redsea_depth_clusters.py). One point
per genome with a resolved depth_m; x = which 70%-identity phaC cluster it
carries, columns ordered by median depth so a depth-segregation pattern, if
real, reads as a trend rather than a shuffled grid. Depth runs top-to-bottom
on y (inverted), matching how depth reads intuitively as a water column.
"""
import csv
import statistics
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

CAT_PALETTE = ['#1B4F9C', '#E0621A', '#1B8A3E', '#C21F6E', '#6A3D9A', '#C2A83E',
               '#1E6E7A', '#9E3B3B', '#5B7C99', '#A8763E', '#3E9E8C', '#556B4F', '#B33951']


def build_depth_cluster_chart(genus, title, out_name, assignments, rows):
    seen = set()
    points = []
    for r in rows:
        if r.get('gtdb_genus', '') != genus or not r.get('depth_m', ''):
            continue
        genome = r['genome']
        if genome in seen:
            continue
        seen.add(genome)
        cid = assignments.get(r.get('target_id', ''))
        if cid:
            points.append((float(r['depth_m']), cid[-12:]))

    by_cluster = {}
    for depth, cid in points:
        by_cluster.setdefault(cid, []).append(depth)

    cluster_order = sorted(by_cluster, key=lambda c: statistics.median(by_cluster[c]))
    x_of = {c: i for i, c in enumerate(cluster_order)}
    colors = {c: CAT_PALETTE[i % len(CAT_PALETTE)] for i, c in enumerate(cluster_order)}

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, ax = plt.subplots(figsize=(9.5, 8), dpi=300)

    for i in range(len(cluster_order)):
        if i % 2 == 0:
            ax.axvspan(i - 0.5, i + 0.5, color='#F2F1EC', zorder=0)

    for c in cluster_order:
        ax.plot([x_of[c], x_of[c]], [min(by_cluster[c]), max(by_cluster[c])],
                 color=colors[c], alpha=0.25, linewidth=6, zorder=1)

    for depth, cid in points:
        x = x_of[cid]
        ax.scatter(x, depth, s=130, color=colors[cid], edgecolor='white', linewidth=0.8, zorder=3, alpha=0.9)

    ax.set_xticks(range(len(cluster_order)))
    ax.set_xticklabels([f'…{c}' for c in cluster_order], family='monospace', fontsize=9, rotation=40, ha='right')
    for tick, c in zip(ax.get_xticklabels(), cluster_order):
        tick.set_color(colors[c])
        tick.set_fontweight('bold')
    ax.set_xlim(-0.7, len(cluster_order) - 0.3)

    ax.set_ylabel('Depth (m)', fontsize=11.5)
    ax.set_xlabel('phaC cluster (70% identity)', fontsize=11.5, labelpad=10)
    ax.invert_yaxis()  # shallow at top, deep at bottom -- reads like a water column
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

    fig.suptitle(title, fontsize=15.5, fontweight='bold', y=0.98)
    fig.text(0.5, 0.925, f'{len(points)} {genus} genomes with resolved depth, {len(cluster_order)} distinct phaC clusters — one point per genome, columns sorted by median depth',
              ha='center', fontsize=9.8, color='#5B6E70')
    fig.text(0.5, 0.01,
              'Faint vertical bars show each cluster\'s observed depth range (min–max). A cluster confined to a narrow band at one\n'
              'end suggests depth-restriction; a cluster spanning the full axis suggests it isn\'t depth-specific.',
              ha='center', va='bottom', fontsize=8.3, color='#5B6E70')

    fig.subplots_adjust(left=0.13, right=0.97, top=0.86, bottom=0.22)

    out_path = OUT / f'{out_name}.png'
    fig.savefig(out_path, dpi=300, facecolor='white')
    print('saved', out_path)
    fig.savefig(OUT / f'{out_name}.pdf', facecolor='white')
    print('saved pdf too')

    for c in cluster_order:
        print(f'{c}: n={len(by_cluster[c])}  depths={sorted(by_cluster[c])}')
