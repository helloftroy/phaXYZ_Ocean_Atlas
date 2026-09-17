"""Companion to figures/phaC_regional_specialists.png -- but the opposite
story. The user asked about two UBA868-dominant PhaC clusters
(...131018176 at 50m, ...095227834 at 90m) and whether there are more.

There are: 23 clusters have UBA868 as their dominant genus. Unlike the
regional-specialist clusters, though, these are NOT geographically
concentrated -- geo_mean_resultant_length for the two named ones is only
~0.30-0.35 (vs >0.85 for every regional specialist), and geo_max_pairwise_km
is ~19,000km, essentially opposite Earth. So a single centroid dot per
cluster (like the regional-specialists figure) would misrepresent them --
the defining fact about these clusters IS that each one is spread across
the whole ocean, not sitting in one place.

So this figure is a small-multiples grid, one global map per cluster,
each showing that cluster's actual per-genome point cloud (not a
centroid) -- restricted to the 10 clusters (of 23 total UBA868-dominant
clusters) with >=50 genomes, to leave out clusters too small/noisy to
show a real spatial pattern.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from phaatlas.pipeline import sequence_clustering as sc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

MIN_GENOMES = 50
ACCENT = '#1E6E7A'

eco_rows = {r['cluster_id']: r for r in csv.DictReader(open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline=''), delimiter='\t')}
selected = [cid for cid, r in eco_rows.items() if r['top_genera'].split(' (')[0] == 'UBA868' and int(r['n_genomes']) >= MIN_GENOMES]
selected.sort(key=lambda c: -int(eco_rows[c]['n_genomes']))
n_total_uba868 = sum(1 for r in eco_rows.values() if r['top_genera'].split(' (')[0] == 'UBA868')

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
points = {c: [] for c in selected}  # (lat, lon) per distinct genome
depths = {c: [] for c in selected}
seen_genomes = {c: set() for c in selected}

with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
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

n_cols, n_rows = 5, 2
fig = plt.figure(figsize=(22, 6.8), dpi=300)

for i, cid in enumerate(selected):
    ax = fig.add_subplot(n_rows, n_cols, i + 1, projection=ccrs.PlateCarree())
    ax.set_global()
    ax.add_feature(cfeature.OCEAN, facecolor='#E4EEEC', zorder=0)
    ax.add_feature(cfeature.LAND, facecolor='#E9E4D6', edgecolor='#B9AF98', linewidth=0.3, zorder=1)
    ax.gridlines(draw_labels=False, linewidth=0.3, color='#C7D0CB', linestyle=(0, (1, 3)), zorder=1)
    for spine in ax.spines.values():
        spine.set_edgecolor('#3A4442')
        spine.set_linewidth(0.6)

    pts = points[cid]
    if pts:
        lats, lons = zip(*pts)
        ax.scatter(lons, lats, transform=ccrs.PlateCarree(), s=7, color=ACCENT, alpha=0.45,
                   linewidths=0, zorder=3, rasterized=True)

    ds = depths[cid]
    mean_depth = f"{sum(ds)/len(ds):.0f} m avg depth" if ds else 'no depth data'
    n_genomes = int(eco_rows[cid]['n_genomes'])
    r_val = float(eco_rows[cid]['geo_mean_resultant_length'])
    ax.set_title(f"…{cid[-12:]}\n{n_genomes:,} genomes  |  {mean_depth}  |  R={r_val:.2f}",
                 fontsize=9.3, fontweight='bold', color='#20302C', pad=5)

fig.suptitle(f'UBA868-Dominant PhaC Clusters Are Cosmopolitan, Not Regional',
             fontsize=18, fontweight='bold', y=0.975)
fig.text(0.5, 0.85,
          f'{len(selected)} of {n_total_uba868} clusters with UBA868 as the dominant host genus, restricted to those with ≥{MIN_GENOMES} genomes — '
          'each panel is one cluster\'s actual genome-level footprint (not a centroid). Contrast with the regional-specialist figure:\n'
          'those clusters had geo_mean_resultant_length > 0.85 (concentrated in one place); every cluster shown here has R < 0.65 — genuinely global, not local.',
          ha='center', va='top', fontsize=10, color='#5B6E70')

fig.subplots_adjust(left=0.02, right=0.98, top=0.70, bottom=0.04, wspace=0.08, hspace=0.55)

out_path = OUT / 'phaC_uba868_global.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phaC_uba868_global.pdf', facecolor='white')
print('saved pdf too')

for cid in selected:
    ds = depths[cid]
    print(f"{cid[-12:]}  n_genomes={eco_rows[cid]['n_genomes']:>5}  n_points={len(points[cid]):>5}  "
          f"mean_depth={sum(ds)/len(ds) if ds else None}  R={eco_rows[cid]['geo_mean_resultant_length']}")
