"""Cluster heatmap (rows = studies, columns = phaC_cluster0.7 size buckets)
answering: does a study's mix of singleton/small/large phaC clusters
track its habitat or geography? Bucket assignment is by a cluster's
GLOBAL size (total distinct genomes carrying it across the whole atlas,
same principle as the rare-cluster tiers in
plot_phac_diversity_habitat_landscape.py) -- so a study's profile shows
whether it mostly rediscovers cosmopolitan clusters or turns up clusters
rare anywhere in the atlas, not just rare within that one study.

Rows: every study with >=50 phaC-positive genomes (100 studies -- the
full ~188-study population was judged too dense to read even with
hierarchical clustering; 50 matches the same MIN_N convention already
used for habitats elsewhere in this project). Hierarchically clustered
(scipy, average-linkage euclidean on the row-normalized proportion
vectors) -- a genuine row dendrogram, not a manually chosen order, so any
grouping by habitat/region in the result is a property of the data, not
the sort. Columns are NOT clustered -- they are a fixed, ordered size
scale (1 / 2-5 / 6-20 / 21-100 / >100), and reordering them would destroy
the one axis that is inherently meaningful.

Values plotted are PROPORTIONS (row-normalized: each study's distinct
clusters / bucket, divided by its own total distinct clusters), not raw
counts -- raw counts are dominated entirely by the handful of
1,000+-genome studies and would make every smaller study's row look
uniformly near-zero, defeating the point of comparing PROFILE SHAPE
across studies of very different sizes. Raw counts are still in the
companion stats TSV for anyone who wants them.

Two renders, per direct request: ALL qualifying studies, and the same
set with "global" (multi-ocean-basin) studies excluded, since a
circumnavigating survey's cluster mix is not really comparable to a
single-site regional study and could dominate/blur any region-level
pattern. "Global" here means longitude span >=180 degrees across the
study's own genomes -- a simple, reproducible proxy (a study confined to
one ocean basin essentially never reaches this; a multi-basin expedition
almost always does). Confirmed against the data: 11/100 studies clear
this bar (TOPC, TPAC, SANC23-1, CAOS20-1, WANG22-1, ...), all large,
well-known multi-region survey-style studies, not an accident of the
threshold.

Two annotation strips between the dendrogram and the heatmap read the
"does it track habitat/region" question directly: each study's DOMINANT
habitat (categorical color, mode of ecosystem_compartment across its own
genomes) and mean |latitude| (sequential color, a simple region/
equator-to-pole proxy -- no ocean-basin classifier exists in this
project to build a real one).

Usage:
    python figures/scripts/plot_phac_cluster_bucket_study_heatmap.py
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

MIN_STUDY_N = 50
GLOBAL_LONSPAN_DEG = 180
BUCKETS = [(1, 1, '1\n(singleton)'), (2, 5, '2–5'), (6, 20, '6–20'),
           (21, 100, '21–100'), (101, float('inf'), '>100')]
EXCLUDE_HABITATS = {  # same list used throughout this project's habitat figures
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}

# ---------------------------------------------------------------------
# 1. genome -> set(cluster_id); cluster -> global size (same construction
#    as plot_phac_diversity_habitat_landscape.py / plot_phac_rank_abundance_rarefaction.py)
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

cluster_genomes = defaultdict(set)
genome_clusters = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        if target_id in bad_targets:
            continue
        for g in target_to_genomes.get(target_id, ()):
            cluster_genomes[cluster_id].add(g)
            genome_clusters[g].add(cluster_id)

cluster_size = {c: len(gs) for c, gs in cluster_genomes.items()}


def bucket_of(size):
    for lo, hi, label in BUCKETS:
        if lo <= size <= hi:
            return label
    raise ValueError(size)


BUCKET_LABELS = [b[2] for b in BUCKETS]
cluster_bucket = {c: bucket_of(s) for c, s in cluster_size.items()}
print(f'{len(cluster_genomes):,} distinct clusters; global bucket sizes: '
      f'{Counter(cluster_bucket.values())}')

# ---------------------------------------------------------------------
# 2. genome -> study/habitat/lat/lon (all from genome_family_matrix.tsv, one file)
# ---------------------------------------------------------------------
genome_study, genome_habitat, genome_lat, genome_lon = {}, {}, {}, {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row['genome']
        if g not in genome_clusters:
            continue
        genome_study[g] = row['study_id']
        h = row['ecosystem_compartment']
        if h not in EXCLUDE_HABITATS:
            genome_habitat[g] = h
        try:
            genome_lat[g] = float(row['latitude_degN'])
            genome_lon[g] = float(row['longitude_degE'])
        except ValueError:
            pass

study_genomes = defaultdict(set)
for g, s in genome_study.items():
    study_genomes[s].add(g)

# ---------------------------------------------------------------------
# 3. per-study bucket profile + habitat/geography annotation
# ---------------------------------------------------------------------
records = []
for s, gs in study_genomes.items():
    if len(gs) < MIN_STUDY_N:
        continue
    clusters = set()
    for g in gs:
        clusters |= genome_clusters[g]
    counts = Counter(cluster_bucket[c] for c in clusters)
    n_clusters = len(clusters)
    lats = [genome_lat[g] for g in gs if g in genome_lat]
    lons = [genome_lon[g] for g in gs if g in genome_lon]
    lonspan = max(lons) - min(lons) if len(lons) > 1 else 0.0
    mean_abs_lat = float(np.mean([abs(x) for x in lats])) if lats else float('nan')
    habs = Counter(genome_habitat.get(g) for g in gs if g in genome_habitat)
    dominant_habitat = habs.most_common(1)[0][0] if habs else 'Mixed/Unknown'
    records.append(dict(
        study=s, n_genomes=len(gs), n_clusters=n_clusters,
        dominant_habitat=dominant_habitat, mean_abs_lat=mean_abs_lat,
        lonspan=lonspan, is_global=lonspan >= GLOBAL_LONSPAN_DEG,
        **{b: counts.get(b, 0) for b in BUCKET_LABELS},
    ))

print(f'{len(records)} studies with >={MIN_STUDY_N} phaC-positive genomes '
      f'({sum(r["is_global"] for r in records)} flagged global, lonspan>={GLOBAL_LONSPAN_DEG}°)')

stats_path = OUT / 'phac_cluster_bucket_study_stats.tsv'
fieldnames = ['study', 'n_genomes', 'n_clusters', 'dominant_habitat', 'mean_abs_lat',
              'lonspan', 'is_global'] + BUCKET_LABELS
with open(stats_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    for r in sorted(records, key=lambda r: -r['n_genomes']):
        writer.writerow({k: (str(r[k]).replace('\n', ' ') if k in BUCKET_LABELS else r[k]) for k in fieldnames})
print(f'wrote {stats_path}')

# ---------------------------------------------------------------------
# categorical habitat palette (dataviz-skill default 8-slot order,
# validated CVD-safe; capped at the top 7 dominant habitats + "Other")
# ---------------------------------------------------------------------
CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
SEQ_BLUES = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
HEATMAP_CMAP = LinearSegmentedColormap.from_list('phac_bucket_seq', SEQ_BLUES, N=256)
LAT_CMAP = LinearSegmentedColormap.from_list('lat_seq', ['#f7f4ee'] + SEQ_BLUES, N=256)
OTHER_COLOR = '#B9B7AE'


def render(records_subset, out_stub, title_suffix):
    if len(records_subset) < 3:
        print(f'skipping {out_stub}: too few studies ({len(records_subset)})')
        return

    habitat_counts = Counter(r['dominant_habitat'] for r in records_subset)
    top_habitats = [h for h, _ in habitat_counts.most_common(len(CATEGORICAL))]
    habitat_color = {h: CATEGORICAL[i] for i, h in enumerate(top_habitats)}

    studies = [r['study'] for r in records_subset]
    matrix_counts = np.array([[r[b] for b in BUCKET_LABELS] for r in records_subset], dtype=float)
    row_sums = matrix_counts.sum(axis=1, keepdims=True)
    matrix_prop = matrix_counts / row_sums

    # hierarchical clustering on row-normalized profiles (rows only -- columns
    # stay in their fixed, meaningful size order)
    Z = linkage(matrix_prop, method='average', metric='euclidean')
    dn = dendrogram(Z, no_plot=True)
    order = dn['leaves']

    matrix_prop = matrix_prop[order]
    studies = [studies[i] for i in order]
    recs_ord = [records_subset[i] for i in order]

    n = len(studies)
    fig_h = max(9, 0.205 * n + 3.2)
    top_margin_in = 1.55   # fixed inches reserved for title+subtitle+rotated strip headers, independent of fig_h
    bottom_margin_in = 0.45
    fig = plt.figure(figsize=(17, fig_h), dpi=300)
    gs = fig.add_gridspec(
        1, 5, width_ratios=[1.1, 0.22, 0.22, 3.6, 0.85],
        wspace=0.04, left=0.035, right=0.80,
        top=1 - top_margin_in / fig_h, bottom=bottom_margin_in / fig_h)
    ax_dendro = fig.add_subplot(gs[0, 0])
    ax_hab = fig.add_subplot(gs[0, 1])
    ax_lat = fig.add_subplot(gs[0, 2])
    ax_heat = fig.add_subplot(gs[0, 3])

    dendrogram(Z, ax=ax_dendro, orientation='left', no_labels=True,
               color_threshold=0, above_threshold_color='#9A9890', link_color_func=lambda k: '#9A9890')
    ax_dendro.invert_yaxis()
    for spine in ax_dendro.spines.values():
        spine.set_visible(False)
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])

    hab_colors = np.array([[matplotlib.colors.to_rgb(habitat_color.get(r['dominant_habitat'], OTHER_COLOR))]
                            for r in recs_ord])
    ax_hab.imshow(hab_colors, aspect='auto', extent=[0, 1, n, 0])
    ax_hab.set_xticks([])
    ax_hab.set_yticks([])
    ax_hab.set_title('Habitat', fontsize=8.5, rotation=90, ha='left', va='bottom', x=0.05, y=1.0)
    for spine in ax_hab.spines.values():
        spine.set_visible(False)

    lat_vals = np.array([r['mean_abs_lat'] for r in recs_ord])
    lat_norm = (lat_vals - np.nanmin(lat_vals)) / (np.nanmax(lat_vals) - np.nanmin(lat_vals) + 1e-9)
    lat_colors = np.array([LAT_CMAP(v)[:3] if np.isfinite(v) else (1, 1, 1) for v in lat_norm]).reshape(n, 1, 3)
    ax_lat.imshow(lat_colors, aspect='auto', extent=[0, 1, n, 0])
    ax_lat.set_xticks([])
    ax_lat.set_yticks([])
    ax_lat.set_title('|Lat|', fontsize=8.5, rotation=90, ha='left', va='bottom', x=0.05, y=1.0)
    for spine in ax_lat.spines.values():
        spine.set_visible(False)

    im = ax_heat.imshow(matrix_prop, aspect='auto', cmap=HEATMAP_CMAP, vmin=0, vmax=1,
                         extent=[0, len(BUCKET_LABELS), n, 0])
    ax_heat.set_xticks(np.arange(len(BUCKET_LABELS)) + 0.5)
    ax_heat.set_xticklabels(BUCKET_LABELS, fontsize=9.5)
    ax_heat.xaxis.set_ticks_position('top')
    ax_heat.set_yticks(np.arange(n) + 0.5)
    label_fontsize = 7.6 if n <= 80 else max(4.8, 7.6 - 0.022 * (n - 80))
    ylabels = [f"{r['study']}  (n={r['n_genomes']:,}{', global' if r['is_global'] else ''})" for r in recs_ord]
    ax_heat.yaxis.tick_right()
    ax_heat.yaxis.set_label_position('right')
    ax_heat.set_yticklabels(ylabels, fontsize=label_fontsize)
    ax_heat.tick_params(left=False, top=False, right=False)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)
    ax_heat.set_xlim(0, len(BUCKET_LABELS))
    ax_heat.set_ylim(n, 0)

    top_frac = 1 - top_margin_in / fig_h
    bottom_frac = bottom_margin_in / fig_h
    cax = fig.add_axes([0.845, bottom_frac + 0.32 * (top_frac - bottom_frac),
                         0.012, 0.22 * (top_frac - bottom_frac)])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Share of that study’s\ndistinct clusters', fontsize=8.7, color='#3A3E3B')
    cbar.ax.tick_params(labelsize=8)
    cbar.outline.set_visible(False)

    legend_handles = [Patch(facecolor=habitat_color[h], edgecolor='none', label=h) for h in top_habitats]
    if any(r['dominant_habitat'] not in habitat_color for r in recs_ord):
        legend_handles.append(Patch(facecolor=OTHER_COLOR, edgecolor='none', label='Other'))
    fig.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(0.815, top_frac),
               fontsize=8.3, frameon=False, title='Dominant habitat', title_fontsize=8.8,
               handlelength=1.1, handleheight=1.1)

    fig.suptitle(f'PhaC cluster-size profile by study{title_suffix}', fontsize=17, fontweight='bold',
                 x=0.04, ha='left', y=1 - 0.42 / fig_h)
    fig.text(0.04, 1 - 0.72 / fig_h,
              f'{n} studies, ≥{MIN_STUDY_N} phaC-positive genomes each. Columns = cluster size in genomes, counted globally across the whole atlas, not just within that study.',
              fontsize=10, color='#5B6E70', ha='left', va='top')

    png_path = OUT / f'{out_stub}.png'
    fig.savefig(png_path, dpi=300, facecolor='white')
    fig.savefig(OUT / f'{out_stub}.pdf', facecolor='white')
    plt.close(fig)
    print(f'saved {png_path} ({n} studies)')


render(records, 'phac_cluster_bucket_study_heatmap_all', ' (all studies)')
non_global = [r for r in records if not r['is_global']]
render(non_global, 'phac_cluster_bucket_study_heatmap_no_global', ' (global surveys excluded)')
