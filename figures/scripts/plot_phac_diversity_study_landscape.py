"""Second version of the PhaC diversity accumulation landscape hero figure
(plot_phac_diversity_habitat_landscape.py): same global rarefaction curve
as sampling-effort baseline, but the overlay is every STUDY (not
habitat) with a phaC-positive genome -- all 188, no minimum-size floor.
(plot_phac_cluster_bucket_study_heatmap.py needed a >=50-genome floor to
keep row labels legible; this is a plain scatter, so no such floor is
needed here -- shown per direct follow-up request, "100 because we only
have 100 or as a choice?".) Colored by each study's DOMINANT habitat
(categorical), not the rare-cluster-enrichment diverging score, per
direct request that the diverging color scheme read as too confusing.
Categorical palette and per-study dominant-habitat logic are the same
method as the study heatmap script, so the color a given habitat gets is
determined the same way in both (most-common-dominant-habitat-first) --
but since this figure's population is all 188 studies rather than the
heatmap's 100, the resulting top-8-habitat ranking (and therefore which
habitats get a named color vs. fold into "Other") is not guaranteed
identical between the two figures.

No permutation-based expected-richness/rare-enrichment stats are
computed here -- this figure is purely observed (n_genomes, n_clusters)
per study against the same global expectation curve; the statistical
habitat-vs-expectation comparison already lives in the habitat-level
hero figure and its stats table.

Usage:
    python figures/scripts/plot_phac_diversity_study_landscape.py
"""
import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

RANDOM_SEED = 42
N_GLOBAL_PERM = 300
N_CURVE_POINTS = 200
MIN_STUDY_N = 2   # excludes 5 single-genome "studies" (3 land exactly on (1,1), overlapping into
                  # what read as one confusing point; 2 more at (1,2)) -- a single genome cannot
                  # really represent a study's diversity, so these are dropped rather than shown

EXCLUDE_HABITATS = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
# same fixed categorical order as plot_phac_cluster_bucket_study_heatmap.py, so
# a habitat's color matches across both figures
CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
OTHER_COLOR = '#B9B7AE'

rng = np.random.default_rng(RANDOM_SEED)

# ---------------------------------------------------------------------
# 1. genome -> set(cluster_id)  (identical construction to the habitat
#    landscape figure and the study heatmap)
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

genome_clusters = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        if target_id in bad_targets:
            continue
        for g in target_to_genomes.get(target_id, ()):
            genome_clusters[g].add(cluster_id)

all_genomes = np.array(sorted(genome_clusters))
n_genomes_total = len(all_genomes)
n_clusters_total = len(set().union(*genome_clusters.values()))
print(f'{n_genomes_total:,} phaC-positive genomes, {n_clusters_total:,} distinct clusters')
genome_cluster_list = [genome_clusters[g] for g in all_genomes]

# ---------------------------------------------------------------------
# 2. genome -> study/habitat, study -> dominant habitat
# ---------------------------------------------------------------------
genome_study, genome_habitat = {}, {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row['genome']
        if g not in genome_clusters:
            continue
        genome_study[g] = row['study_id']
        h = row['ecosystem_compartment']
        if h not in EXCLUDE_HABITATS:
            genome_habitat[g] = h

study_genome_idx = defaultdict(list)
for i, g in enumerate(all_genomes):
    s = genome_study.get(g)
    if s:
        study_genome_idx[s].append(i)

study_records = []
for s, idx in study_genome_idx.items():
    if len(idx) < MIN_STUDY_N:
        continue
    clusters = set()
    for i in idx:
        clusters |= genome_cluster_list[i]
    habs = Counter(genome_habitat.get(all_genomes[i]) for i in idx if all_genomes[i] in genome_habitat)
    dominant = habs.most_common(1)[0][0] if habs else 'Mixed/Unknown'
    study_records.append(dict(study=s, n_genomes=len(idx), n_clusters=len(clusters), dominant_habitat=dominant))

print(f'{len(study_records)} studies with >={MIN_STUDY_N} phaC-positive genomes')
habitat_counts = Counter(r['dominant_habitat'] for r in study_records)
top_habitats = [h for h, _ in habitat_counts.most_common(len(CATEGORICAL))]
habitat_color = {h: CATEGORICAL[i] for i, h in enumerate(top_habitats)}
print('dominant-habitat counts among these studies:', habitat_counts.most_common())

# ---------------------------------------------------------------------
# 3. global rarefaction curve (identical method/seed to the habitat landscape figure)
# ---------------------------------------------------------------------
def downsample_indices(n, n_points):
    if n <= n_points:
        return np.arange(1, n + 1)
    return np.unique(np.round(np.linspace(1, n, n_points)).astype(int))


xs_global = downsample_indices(n_genomes_total, N_CURVE_POINTS)
curves = np.zeros((N_GLOBAL_PERM, len(xs_global)), dtype=float)
order_template = np.arange(n_genomes_total)
for p in range(N_GLOBAL_PERM):
    order = rng.permutation(order_template)
    seen = set()
    xi = 0
    for i, gi in enumerate(order, start=1):
        seen.update(genome_cluster_list[gi])
        if xi < len(xs_global) and i == xs_global[xi]:
            curves[p, xi] = len(seen)
            xi += 1
print(f'global rarefaction: {N_GLOBAL_PERM} permutations done')
mean_curve = curves.mean(axis=0)
lo_curve = np.percentile(curves, 2.5, axis=0)
hi_curve = np.percentile(curves, 97.5, axis=0)

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5})
CURVE_COLOR = '#3A4A46'
CURVE_BAND = '#E4E8E5'
TEXT_MUTED = '#5B6E70'
EDGE_WHITE = '#FCFCFB'

fig, ax = plt.subplots(figsize=(13.5, 9.8), dpi=300)

ax.fill_between(xs_global, lo_curve, hi_curve, color=CURVE_BAND, zorder=1, linewidth=0,
                 label='95% permutation interval')
ax.plot(xs_global, mean_curve, color=CURVE_COLOR, linewidth=2.2, zorder=2,
        label='Global expectation (mean of\nrandom genome subsamples)')

for h in top_habitats:
    pts = [r for r in study_records if r['dominant_habitat'] == h]
    ax.scatter([r['n_genomes'] for r in pts], [r['n_clusters'] for r in pts],
               s=95, color=habitat_color[h], edgecolors=EDGE_WHITE, linewidths=1.1,
               zorder=5, label=None)
other_pts = [r for r in study_records if r['dominant_habitat'] not in habitat_color]
if other_pts:
    ax.scatter([r['n_genomes'] for r in other_pts], [r['n_clusters'] for r in other_pts],
               s=95, color=OTHER_COLOR, edgecolors=EDGE_WHITE, linewidths=1.1, zorder=5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('PhaC-positive genomes sampled (log scale)', fontsize=12)
ax.set_ylabel('Distinct phaC_cluster0.7 clusters recovered (log scale)', fontsize=12)

# tight axis limits: the plot's origin sits flush at the data's own floor
# (no default matplotlib log-margin padding below/left of it), with a
# little headroom only above/right of the data for labels/legend to sit in
x_floor = min(xs_global[0], min(r['n_genomes'] for r in study_records))
y_floor = min(mean_curve[0], min(r['n_clusters'] for r in study_records))
x_ceil = max(xs_global[-1], max(r['n_genomes'] for r in study_records))
y_ceil = max(hi_curve[-1], max(r['n_clusters'] for r in study_records))
ax.set_xlim(x_floor, x_ceil * 1.15)
ax.set_ylim(y_floor, y_ceil * 1.25)
ax.grid(True, which='major', color='#EBEEEC', linewidth=0.6, zorder=0)
ax.grid(True, which='minor', color='#F4F6F4', linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('#C3C2B7')
ax.spines['bottom'].set_color('#C3C2B7')

habitat_legend = [Patch(facecolor=habitat_color[h], edgecolor='none', label=h) for h in top_habitats]
if other_pts:
    habitat_legend.append(Patch(facecolor=OTHER_COLOR, edgecolor='none', label='Other'))
leg1 = ax.legend(handles=habitat_legend, loc='lower right', fontsize=9, frameon=False,
                  title='Dominant habitat', title_fontsize=9.7, ncol=1, handlelength=1.1)
ax.add_artist(leg1)
# curve/band legend separately (upper-left, away from the point cloud)
band_patch = Patch(facecolor=CURVE_BAND, edgecolor='none', label='95% permutation interval')
line_handle = plt.Line2D([0], [0], color=CURVE_COLOR, linewidth=2.2,
                          label='Global expectation (mean of\nrandom genome subsamples)')
ax.legend(handles=[line_handle, band_patch], loc='upper left', fontsize=9, frameon=False, handlelength=2.2)

stats_text = (
    f'{n_clusters_total:,} clusters from {n_genomes_total:,} phaC-positive genomes\n'
    f'all {len(study_records)} studies with a phaC-positive genome\n'
    f'colored by that study’s own dominant habitat'
)
ax.text(0.015, 0.72, stats_text, transform=ax.transAxes, ha='left', va='top', fontsize=9.3,
        color='#3A4A46', linespacing=1.55,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F7F8F6', edgecolor='#D8DCD9', linewidth=0.8))

fig.suptitle('The PhaC diversity accumulation landscape, by study', fontsize=19, fontweight='bold',
              y=0.985, x=0.065, ha='left')
fig.text(0.065, 0.945,
         'Every study with a phaC-positive genome, plotted against the same global sampling-effort curve.\n'
         'Color identifies each study’s dominant habitat -- position above/below the curve still reads richness vs. expectation.',
         ha='left', va='top', fontsize=10.5, color=TEXT_MUTED, linespacing=1.5)

fig.tight_layout(rect=[0.005, 0.01, 0.99, 0.90])

out_path = OUT / 'phac_diversity_study_landscape.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_diversity_study_landscape.pdf', facecolor='white')
