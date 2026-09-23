"""Single-map overlay of a genus's largest globally-distributed phaC
clusters -- generalized from the earlier UBA868-specific overlay after
that genus turned out to be a middling choice post-audit (see
plot_genus_global_clusters.py and figures/PHA_CLEAN_RESULTS.md section
5.2 for why). Same per-genome point-cloud approach as the small-multiples
version (these clusters are cosmopolitan, not regional, so a single
centroid dot per cluster would misrepresent them), but overlaid on one
map with a per-cluster legend showing average depth (mean depth_m across
distinct genomes with resolved depth, missing dropped), matching the
convention set in plot_regional_specialists.py.

QC: excludes rows whose best_query is one of 67 confirmed-wrong-gene phaC
reference proteins -- see _phac_qc.py. phaC_cluster0.7_cluster_ecology.tsv
was regenerated with the same filter applied in pipeline/cluster_ecology.py.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from phaatlas.pipeline import sequence_clustering as sc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
OUT = Path(__file__).resolve().parent.parent

MIN_GENOMES = 200
CAT_PALETTE = ['#1B4F9C', '#E0621A', '#1B8A3E', '#C21F6E', '#6A3D9A', '#C2A83E', '#1E6E7A']
MARKERS = ['o', '^', 's', 'D', 'v', 'P', 'X']

GENUS = sys.argv[1] if len(sys.argv) > 1 else 'Sulfitobacter'
# same fix as plot_genus_global_clusters.py -- don't mangle GTDB polyphyletic
# suffixes like "Pseudomonas_E" into a nonsense word
OUT_STEM = f'phaC_{GENUS.lower()}_overlay'

eco_rows = {r['cluster_id']: r for r in csv.DictReader(open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline=''), delimiter='\t')}
n_total_genus = sum(1 for r in eco_rows.values() if r['top_genera'].split(' (')[0] == GENUS)
selected = [cid for cid, r in eco_rows.items() if r['top_genera'].split(' (')[0] == GENUS and int(r['n_genomes']) > MIN_GENOMES]
selected.sort(key=lambda c: -int(eco_rows[c]['n_genomes']))
colors = {cid: CAT_PALETTE[i % len(CAT_PALETTE)] for i, cid in enumerate(selected)}
markers = {cid: MARKERS[i % len(MARKERS)] for i, cid in enumerate(selected)}

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
points = {c: [] for c in selected}
depths = {c: [] for c in selected}
seen_genomes = {c: set() for c in selected}

with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if _phac_qc.is_bad(row.get('best_query', '')):
            continue
        cid = assignments.get(row.get('target_id', ''))
        if cid not in points:
            continue
        genome = row.get('genome', '')
        if not genome or genome in seen_genomes[cid]:
            continue
        seen_genomes[cid].add(genome)
        lat, lon = row.get('latitude_degN', ''), row.get('longitude_degE', '')
        if lat and lon:
            points[cid].append((float(lat), float(lon)))
        d = row.get('depth_m', '')
        if d:
            depths[cid].append(float(d))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})

fig = plt.figure(figsize=(19, 9.2), dpi=300)
MAP_LEFT = 0.36
ax = fig.add_axes([MAP_LEFT, 0.09, 0.99 - MAP_LEFT, 0.76], projection=ccrs.PlateCarree())
ax.set_global()
ax.add_feature(cfeature.OCEAN, facecolor='#E4EEEC', zorder=0)
ax.add_feature(cfeature.LAND, facecolor='#E9E4D6', edgecolor='#B9AF98', linewidth=0.5, zorder=1)
ax.gridlines(draw_labels=False, linewidth=0.4, color='#C7D0CB', linestyle=(0, (1, 3)), zorder=1)
for spine in ax.spines.values():
    spine.set_edgecolor('#3A4442')
    spine.set_linewidth(0.9)

for cid in sorted(selected, key=lambda c: -int(eco_rows[c]['n_genomes'])):
    pts = points[cid]
    if not pts:
        continue
    lats, lons = zip(*pts)
    ax.scatter(lons, lats, transform=ccrs.PlateCarree(), s=26, color=colors[cid], alpha=0.8,
               marker=markers[cid], edgecolor='white', linewidths=0.3,
               zorder=2 + selected.index(cid), rasterized=True)

fig.suptitle(f'{GENUS}: Do Its Largest PhaC Clusters Share a Range?', fontsize=18, fontweight='bold',
             x=(MAP_LEFT + 0.99) / 2, y=0.975)
fig.text((MAP_LEFT + 0.99) / 2, 0.905,
          f'The {len(selected)} {GENUS}-dominant clusters with >{MIN_GENOMES} genomes (of {n_total_genus} total), overlaid — one point per distinct genome, colored by cluster',
          ha='center', fontsize=10.5, color='#5B6E70')

legend_lines = []
for cid in sorted(selected, key=lambda c: (depths[c] == [], -(sum(depths[c]) / len(depths[c])) if depths[c] else 0)):
    ds = depths[cid]
    depth_str = f"{sum(ds)/len(ds):.0f} m (n={len(ds)})" if ds else 'no depth data'
    legend_lines.append((colors[cid], markers[cid], cid[-12:], int(eco_rows[cid]['n_genomes']), depth_str))

legend_x, legend_y0, dy = 0.025, 0.12, 0.075
fig.text(legend_x, legend_y0 + dy * len(legend_lines) + 0.035, 'Cluster', fontsize=11.5, fontweight='bold', color='#20302C', va='bottom')
fig.text(legend_x + 0.155, legend_y0 + dy * len(legend_lines) + 0.035, 'Avg. depth', fontsize=11.5, fontweight='bold', color='#20302C', va='bottom')
legend_marker_ax = fig.add_axes([0, 0, 1, 1], zorder=20)
legend_marker_ax.axis('off')
legend_marker_ax.set_xlim(0, 1)
legend_marker_ax.set_ylim(0, 1)
for i, (color, marker, short_id, n_genomes, depth_str) in enumerate(legend_lines):
    y = legend_y0 + dy * (len(legend_lines) - 1 - i)
    legend_marker_ax.scatter([legend_x + 0.01], [y], s=170, color=color, marker=marker,
                              edgecolor='white', linewidths=0.6, transform=fig.transFigure, clip_on=False)
    fig.text(legend_x + 0.03, y + 0.012, f'…{short_id}', fontsize=10.5, fontweight='bold', va='center', family='monospace')
    fig.text(legend_x + 0.03, y - 0.012, f'{n_genomes:,} genomes', fontsize=9, color='#5B6E70', va='center')
    fig.text(legend_x + 0.155, y, depth_str, fontsize=10, color='#5B6E70', va='center')

footnote = (
    "Average depth = mean depth_m across distinct genomes with resolved depth in that cluster (missing values dropped, not imputed as\n"
    "zero). Circles plotted largest-cluster-first so smaller clusters' points stay visible on top; overlap between clusters is real overlap."
)
fig.text((MAP_LEFT + 0.99) / 2, 0.015, footnote, ha='center', va='bottom', fontsize=8.4, color='#5B6E70')

out_path = OUT / f'{OUT_STEM}.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / f'{OUT_STEM}.pdf', facecolor='white')
print('saved pdf too')

for cid in selected:
    ds = depths[cid]
    print(f"{cid[-12:]}  n_genomes={eco_rows[cid]['n_genomes']:>5}  mean_depth={sum(ds)/len(ds) if ds else None}")
