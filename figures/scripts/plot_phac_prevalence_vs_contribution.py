"""Prevalence vs. contribution, one bubble per phylum.
x = % of that phylum's sampled genomes carrying phaC (prevalence)
y = that phylum's share of all phaC-positive genomes (log scale)
size = distinct phaC_cluster0.7 clusters found in the phylum

Separates "common within the lineage" from "supplies most of the dataset":
Pseudomonadota is high on y, Myxococcota is far right on x.

Denominators as in plot_phac_prevalence_by_phylum.py (all 274,282 OMDB
genomes, omdb_all_genome_taxonomy.tsv). Phyla with >=500 sampled genomes
and at least one phaC-positive genome.

Usage:
    python figures/scripts/plot_phac_prevalence_vs_contribution.py
"""
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
MIN_GENOMES = 500
# Bubbles take the phylum's own colour from _phylum_colors, the same one the
# radial tree and the synthase-class bars use, because the three are assembled
# into one taxonomy figure and a phylum has to be identifiable across all of
# them. The previous teal/orange split encoded only "above or below the
# dataset-wide prevalence", which the dashed reference line already says, and
# spent the figure's one free channel saying it twice.
import _phylum_colors  # noqa: E402
TEXT_DARK, TEXT_MUTED, GRID = '#20302C', '#5B6E70', '#E9ECEA'

phac_pos = set()
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if int(row['n_phaC']) > 0:
            phac_pos.add(row['genome'])

phylum_of, total, pos = {}, Counter(), Counter()
n_all = 0
with open(FA / 'omdb_all_genome_taxonomy.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n_all += 1
        p = row['gtdb_phylum']
        if not p or p.lower().startswith('unknown'):
            continue
        phylum_of[row['genome']] = p
        total[p] += 1
        if row['genome'] in phac_pos:
            pos[p] += 1
overall = 100 * len(phac_pos) / n_all

# distinct clusters per phylum (QC-passing targets, canonical phaC-positive genomes)
bad = _phac_qc.load_bad_targets()
target_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for c, t in csv.reader(f, delimiter='\t'):
        if t not in bad:
            target_cluster[t] = c
clusters = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t'); next(r)
    for t, g in r:
        if t in target_cluster and g in phac_pos and g in phylum_of:
            clusters[phylum_of[g]].add(target_cluster[t])

rows = []
for p, n in total.items():
    if n < MIN_GENOMES or pos[p] == 0:
        continue
    rows.append(dict(phylum=p, n_genomes=n, n_phac=pos[p], prevalence=100 * pos[p] / n,
                     share=100 * pos[p] / len(phac_pos), n_clusters=len(clusters[p])))
rows.sort(key=lambda r: -r['share'])
with open(OUT / 'phac_prevalence_vs_contribution.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f'{v:.3f}' if isinstance(v, float) else v) for k, v in r.items()})
        print(f"{r['phylum']:20s} prevalence {r['prevalence']:5.1f}%  share {r['share']:6.2f}%  clusters {r['n_clusters']:5,}")

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=300)
max_c = max(r['n_clusters'] for r in rows)


def size_of(c):
    return 28 + 1700 * c / max_c


ax.axvline(overall, color=TEXT_DARK, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
for r in sorted(rows, key=lambda r: -r['n_clusters']):
    ax.scatter(r['prevalence'], r['share'], s=size_of(r['n_clusters']),
               color=_phylum_colors.color_for(r['phylum']), alpha=0.88,
               edgecolor='white', linewidth=1.0, zorder=3)

LABEL = {'Pseudomonadota': (0, -30, 'center'), 'Myxococcota': (0, 12, 'center'), 'Thermoproteota': (12, 6, 'left'),
         'SAR324': (10, -10, 'left'), 'Desulfobacterota': (10, -14, 'left'), 'Actinomycetota': (10, 4, 'left'),
         'Chloroflexota': (11, 0, 'left'), 'Acidobacteriota': (10, -2, 'left'), 'Bacillota': (9, -4, 'left'),
         'Bacteroidota': (9, 9, 'left'), 'Cyanobacteriota': (10, -9, 'left')}
for r in rows:
    if r['phylum'] in LABEL:
        dx, dy, ha = LABEL[r['phylum']]
        ax.annotate(r['phylum'], (r['prevalence'], r['share']), xytext=(dx, dy), textcoords='offset points',
                    ha=ha, va='center', fontsize=8.8, color=TEXT_DARK, zorder=5,
                    path_effects=[pe.withStroke(linewidth=2.4, foreground='white')])

ax.set_yscale('log')
ax.set_ylim(0.002, 400)
ax.set_xlim(-2, 44)
ax.set_yticks([0.01, 0.1, 1, 10, 100])
ax.set_yticklabels(['0.01%', '0.1%', '1%', '10%', '100%'])
ax.set_xlabel('Prevalence: genomes in the phylum carrying phaC (%)', fontsize=10.5, color=TEXT_DARK)
ax.set_ylabel('Contribution: share of all phaC-positive genomes (log)', fontsize=10.5, color=TEXT_DARK)
ax.grid(True, which='major', color=GRID, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color('#C3C2B7')
ax.tick_params(colors=TEXT_MUTED, length=3)
ax.text(overall + 0.5, 0.004, f'all genomes: {overall:.1f}%', fontsize=8.5, color=TEXT_DARK, va='center')

size_handles = [plt.scatter([], [], s=size_of(c), color='#B9C0BC', edgecolor='white', label=f'{c:,}') for c in (50, 500, 2500)]
ax.legend(handles=size_handles, title='phaC clusters', title_fontsize=8.6, fontsize=8.4, frameon=False,
          loc='lower right', ncol=3, columnspacing=2.6, borderpad=1.4, handletextpad=1.6)
fig.tight_layout()
fig.savefig(OUT / 'phac_prevalence_vs_contribution.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'phac_prevalence_vs_contribution.pdf', facecolor='white')
print('saved', OUT / 'phac_prevalence_vs_contribution.png')
