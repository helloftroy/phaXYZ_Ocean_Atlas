"""Small-multiples world-map grid of a genus's globally-distributed phaC
clusters -- generalized from the earlier UBA868-specific figure after
that genus turned out to be a middling choice post-audit (only 6
qualifying clusters, 711 genomes) once the two clusters that originally
prompted picking it were found to be near-total contamination (see
figures/PHA_CLEAN_RESULTS.md section 5.1). Re-scanned every genus in the
corrected phaC_cluster0.7_cluster_ecology.tsv for ones with multiple
large (>=30 genomes), genuinely globally-spread (geo_mean_resultant_length
< 0.4, max pairwise span > 10,000km) clusters -- three real, well-known
marine genera came out on top: Marinobacter, Sulfitobacter, Pseudomonas_E.
This script builds the same figure for any one of them, run three times.

QC: excludes rows whose best_query is one of 67 confirmed-wrong-gene phaC
reference proteins -- see _phac_qc.py. phaC_cluster0.7_cluster_ecology.tsv
(cluster selection/counts) was regenerated with the same filter applied in
pipeline/cluster_ecology.py.
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

MIN_GENOMES = 30
ACCENT = '#1E6E7A'

GENUS = sys.argv[1] if len(sys.argv) > 1 else 'Marinobacter'
# NOT .replace("_", "") -- that mangled GTDB's polyphyletic-split suffix (e.g.
# "Pseudomonas_E") into a nonsense word ("pseudomonase") and silently wrote a
# second, wrongly-named file instead of updating the one this project's docs
# reference (phaC_pseudomonas_e_global.png). Found while rerunning this script
# post the 2026-09-22 phaC fix.
OUT_STEM = f'phaC_{GENUS.lower()}_global'

MAX_R = 0.4  # "genuinely global, not regional" threshold -- must match the footnote text below

def is_global(r):
    if not r['geo_mean_resultant_length']:
        return False
    return float(r['geo_mean_resultant_length']) < MAX_R

eco_rows = {r['cluster_id']: r for r in csv.DictReader(open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline=''), delimiter='\t')}
genus_rows = {cid: r for cid, r in eco_rows.items() if r['top_genera'].split(' (')[0] == GENUS}
n_total_genus = len(genus_rows)
selected = [cid for cid, r in genus_rows.items() if int(r['n_genomes']) >= MIN_GENOMES and is_global(r)]
selected.sort(key=lambda c: -int(eco_rows[c]['n_genomes']))
n_qualifying_by_size = sum(1 for r in genus_rows.values() if int(r['n_genomes']) >= MIN_GENOMES)

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
points = {c: [] for c in selected}  # (lat, lon) per distinct genome
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

n = len(selected)
n_cols = min(5, n) if n else 1
n_rows = max(1, -(-n // n_cols))
fig = plt.figure(figsize=(4.4 * n_cols, 3.6 * n_rows + 1.6), dpi=300)

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

fig.suptitle(f'{GENUS}-Dominant PhaC Clusters Are Cosmopolitan, Not Regional',
             fontsize=18, fontweight='bold', y=0.975)
fig.text(0.5, 0.90,
          f'{len(selected)} of {n_total_genus} clusters with {GENUS} as the dominant host genus qualify: ≥{MIN_GENOMES} genomes AND geo_mean_resultant_length < {MAX_R} '
          f'(of {n_qualifying_by_size} large enough by genome count alone, several were more regionally concentrated and are excluded here). '
          'Each panel is one cluster\'s actual genome-level footprint (not a centroid).\n'
          'Contrast with the regional-specialist figure: those clusters had geo_mean_resultant_length > 0.85 (concentrated in one place); every cluster shown here has R < 0.4 — genuinely global, not local.',
          ha='center', va='top', fontsize=9.5, color='#5B6E70')

fig.subplots_adjust(left=0.02, right=0.98, top=1 - 1.6 / (3.6 * n_rows + 1.6), bottom=0.04, wspace=0.08, hspace=0.55)

out_path = OUT / f'{OUT_STEM}.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / f'{OUT_STEM}.pdf', facecolor='white')
print('saved pdf too')

for cid in selected:
    ds = depths[cid]
    print(f"{cid[-12:]}  n_genomes={eco_rows[cid]['n_genomes']:>5}  n_points={len(points[cid]):>5}  "
          f"mean_depth={sum(ds)/len(ds) if ds else None}  R={eco_rows[cid]['geo_mean_resultant_length']}")
