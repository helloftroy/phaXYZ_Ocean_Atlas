"""Two 100%-stacked bars communicating the phaC_cluster0.7 long tail more
directly than the rank-abundance panel (plot_phac_rank_abundance_rarefaction.py
Panel A) does: the same 5 occupancy buckets as
plot_phac_cluster_bucket_study_heatmap.py (1 / 2-5 / 6-20 / 21-100 / >100
genomes per cluster, bucketed by a cluster's GLOBAL size), but shown as
the share of ALL distinct clusters in each bucket (top bar) against the
share of all genome-cluster MEMBERSHIPS -- "occupancy," i.e. summing
cluster size itself, not counting clusters -- in each bucket (bottom
bar). The two bars are deliberately near-mirror images: this is the
point. Most cluster DIVERSITY sits in the low-occupancy buckets; most
actual phaC gene copies sit in a small number of high-occupancy
clusters. Neither bar alone tells the whole story.

"Occupancy" here means the sum of cluster sizes (genome-cluster
membership edges), not distinct genome count -- a genome with 2 phaC
paralogs in the same cluster contributes 2 to that cluster's size, same
convention every other cluster-size figure in this project uses. Labeled
explicitly as "phaC gene copies" in the figure to avoid reading as "share
of genomes," which it is not.

Usage:
    python figures/scripts/plot_phac_cluster_occupancy_spectrum.py
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

BUCKETS = [(1, 1, '1 genome\n(singleton)'), (2, 5, '2–5'), (6, 20, '6–20'),
           (21, 100, '21–100'), (101, float('inf'), '>100')]
# 5-step ordinal ramp from the same sequential-blue family used throughout
# this session's figures (heatmap, landscape plots), light -> dark
BUCKET_COLORS = ['#9ec5f4', '#6da7ec', '#2a78d6', '#1c5cab', '#0d366b']

# ---------------------------------------------------------------------
# 1. cluster -> global size (identical construction to every other
#    cluster-size figure in this project)
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
target_to_genomes = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        if target_id in bad_targets:
            continue
        target_to_genomes[target_id].add(genome)

# Genome pool is the full protein-to-genome membership list above, QC-filtered --
# 38,744 phaC-positive genomes. An earlier version restricted this to the 31,464
# genomes with n_phaC > 0 in genome_family_matrix.tsv, on the understanding that
# the larger number was an unexplained join artifact. It was the opposite: the
# matrix comes from a metadata table capped at 5 genomes per identical protein,
# which hid 7,298 real phaC carriers (PHA_CLEAN_RESULTS.md section 2). The
# uncapped list is the correct pool, so no second filter is applied here.
cluster_genomes = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        if target_id in bad_targets:
            continue
        for g in target_to_genomes.get(target_id, ()):
            cluster_genomes[cluster_id].add(g)
cluster_genomes = {c: g for c, g in cluster_genomes.items() if g}
n_genomes = len(set().union(*cluster_genomes.values()))

cluster_size = {c: len(gs) for c, gs in cluster_genomes.items()}
n_clusters_total = len(cluster_size)
total_occupancy = sum(cluster_size.values())
print(f'{n_clusters_total:,} distinct clusters, {total_occupancy:,} total genome-cluster memberships')


def bucket_of(size):
    for lo, hi, label in BUCKETS:
        if lo <= size <= hi:
            return label
    raise ValueError(size)


BUCKET_LABELS = [b[2] for b in BUCKETS]
cluster_count = Counter()
occupancy = Counter()
for c, s in cluster_size.items():
    b = bucket_of(s)
    cluster_count[b] += 1
    occupancy[b] += s

cluster_pct = [100 * cluster_count[b] / n_clusters_total for b in BUCKET_LABELS]
occupancy_pct = [100 * occupancy[b] / total_occupancy for b in BUCKET_LABELS]
for label, cp, op in zip(BUCKET_LABELS, cluster_pct, occupancy_pct):
    print(f"{label.replace(chr(10),' '):18s} clusters={cp:5.1f}%   genome-cluster pairs={op:5.1f}%")

stats_path = OUT / 'phac_cluster_occupancy_spectrum_stats.tsv'
with open(stats_path, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['bucket', 'n_clusters', 'pct_of_clusters', 'n_genome_cluster_memberships', 'pct_of_occupancy'])
    for label in BUCKET_LABELS:
        b = label.replace('\n', ' ')
        writer.writerow([b, cluster_count[label], round(100 * cluster_count[label] / n_clusters_total, 2),
                          occupancy[label], round(100 * occupancy[label] / total_occupancy, 2)])
print(f'wrote {stats_path}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12})
TEXT_MUTED = '#5B6E70'
TEXT_DARK = '#1E2422'

fig, axes = plt.subplots(2, 1, figsize=(12, 4.6), dpi=300)
bar_rows = [
    (axes[0], cluster_pct, 'Share of distinct\nphaC clusters', f'{n_clusters_total:,} clusters'),
    (axes[1], occupancy_pct, 'Share of genome\u2013cluster pairs\n(genome counted once per cluster)', f'{total_occupancy:,} pairs\nin {n_genomes:,} genomes'),
]

for ax, values, row_label, n_label in bar_rows:
    left = 0
    n_small = 0
    for val, color, blabel in zip(values, BUCKET_COLORS, BUCKET_LABELS):
        ax.barh([0], [val], left=left, color=color, height=0.62, edgecolor='#FCFCFB', linewidth=1.5)
        # label inside the segment if it's wide enough, else above with a leader
        text_color = '#FFFFFF' if color in ('#2a78d6', '#1c5cab', '#0d366b') else TEXT_DARK
        if val >= 9:
            ax.text(left + val / 2, 0, f'{val:.1f}%', ha='center', va='center',
                     fontsize=11.5, color=text_color, fontweight='medium')
        else:
            ax.annotate(f'{val:.1f}%', xy=(left + val / 2, 0.31), xytext=(left + val / 2, 0.62 + 0.3 * n_small),
                        ha='center', va='bottom', fontsize=9.7, color=TEXT_DARK,
                        arrowprops=dict(arrowstyle='-', color='#9A9890', linewidth=0.8))
        n_small = n_small + 1 if val < 9 else 0
        left += val
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.55, 1.25)
    ax.set_yticks([0])
    ax.set_yticklabels([row_label], fontsize=11, ha='right')
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(100.5, 0, n_label, ha='left', va='center', fontsize=9.5, color=TEXT_MUTED)

legend_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor='none')
                  for c in BUCKET_COLORS]
fig.legend(legend_handles, [b.replace('\n', ' ') for b in BUCKET_LABELS],
           loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.02), frameon=False,
           fontsize=10.3, title='Cluster size (genomes carrying it, counted across the whole atlas)',
           title_fontsize=10, handlelength=1.3, columnspacing=1.6)

# No title or subtitle: this figure goes into a manuscript, where the caption is
# set by the journal and anything written inside the image gets cropped off.
fig.tight_layout(rect=[0.10, 0.06, 0.88, 1.0])

out_path = OUT / 'phac_cluster_occupancy_spectrum.png'
fig.savefig(out_path, dpi=300, facecolor='white', bbox_inches='tight')
print('saved', out_path)
fig.savefig(OUT / 'phac_cluster_occupancy_spectrum.pdf', facecolor='white', bbox_inches='tight')
