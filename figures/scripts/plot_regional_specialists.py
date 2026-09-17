"""Global map of PhaC 70%-identity clusters that are geographically
"regional specialists" -- concentrated in one place/region rather than
cosmopolitan -- rather than the whole 20,211-cluster set, which would be
unplottable as a map.

Selection: from fair_ocean_agent/phaC_cluster0.7_cluster_ecology.tsv,
candidates need geo_mean_resultant_length > 0.85 (a circular-statistics
concentration measure: 1.0 = every genome pulls the same compass direction
from Earth's center, i.e. genuinely one place; 0 = uniformly scattered) and
enough genomes with resolved depth (n_genomes_with_depth >= 8) to make the
"average depth" legend entry meaningful -- most candidates have zero depth
coverage (depth metadata is sparse project-wide), so this filter matters.
The final 12 were then hand-picked from that candidate pool for global
spread and biological interest (surface vs. abyssal, vent-associated
archaea, Arctic vs. tropical), not just top-N by genome count.

Two things this script does NOT take from the ecology TSV at face value:
  - Position: uses the TSV's own geo_centroid_lat/lon (a proper circular
    mean), not a naive arithmetic mean of raw lat/lon -- confirmed live
    that naively averaging longitudes for a cluster with any real spread
    silently produces a wrong centroid (a few of these candidates showed
    >5,000km position errors this way, e.g. across the Pacific vs.
    Atlantic, from longitude wraparound).
  - Average depth: recomputed directly from
    phaC_unique_targets_with_metadata_depth.tsv per distinct genome,
    dropping genomes with no resolved depth, per the user's explicit ask
    ("average after dropping missing depths") -- the ecology TSV's own
    median_depth_m is a median, a different statistic, not reused here.

Overlap handling: any two selected clusters within OVERLAP_KM of each other
(only the Sydney, Australia pair here) are drawn as a small true-position
anchor dot connected by a thin leader line to an offset, larger label
marker -- otherwise two colored dots would sit exactly on top of each
other and be unreadable.
"""
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from phaatlas.pipeline import sequence_clustering as sc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

# hand-picked from the candidate pool (see module docstring) for global
# spread + biological interest; the Sydney pair is deliberately kept in to
# demonstrate the overlap/leader-line handling this figure needed.
SELECTED = [
    '000033545441',  # SHLQ01, Sydney AU -- overlaps with the next one
    '000034200129',  # SHLQ01, Sydney AU
    '000140439260',  # Moritella, NE Pacific, abyssal
    '000067604172',  # Thermococcus/Pyrococcus, S. Atlantic, hydrothermal-vent depth
    '000045952813',  # UBA4427, sub-Arctic (Davis Strait)
    '000005303736',  # UBA4582, N. Atlantic (Norwegian shelf)
    '000232749406',  # Crocosphaera, tropical N. Pacific
    '000012716247',  # JACZQZ01, NW Pacific (Philippine Sea)
    '000226102035',  # Robiginitomaculum_A, mid N. Atlantic
    '000018400463',  # Rs1, N. Pacific (Hawaii), abyssal
    '000001867111',  # BACL27, Baltic Sea
    '000027486053',  # Algiphilus, Southern Ocean (south of Australia)
]

OVERLAP_KM = 300
OFFSET_DEG = 7.0

CAT_PALETTE = ['#1E6E7A', '#C9622D', '#8B5FBF', '#3E8914', '#C2A83E', '#B33951',
               '#4A7FB5', '#D98E04', '#6B4226', '#6FA88A', '#9E3B3B', '#7A7A7A']


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


# ---- load ecology rows for the selected clusters ----
eco_rows = {r['cluster_id'][-12:]: r for r in csv.DictReader(open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline=''), delimiter='\t')}
clusters = []
for short_id in SELECTED:
    r = eco_rows[short_id]
    genus = r['top_genera'].split(' (')[0]
    clusters.append({
        'id': short_id,
        'full_id': r['cluster_id'],
        'lat': float(r['geo_centroid_lat']),
        'lon': float(r['geo_centroid_lon']),
        'n_genomes': int(r['n_genomes']),
        'genus': genus,
    })

# ---- recompute mean depth per cluster, dropping missing, from raw metadata ----
assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
full_ids = {c['full_id'] for c in clusters}
depths = {fid: [] for fid in full_ids}
seen_genomes = {fid: set() for fid in full_ids}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cid = assignments.get(row.get('target_id', ''))
        if cid not in full_ids:
            continue
        genome = row.get('genome', '')
        if not genome or genome in seen_genomes[cid]:
            continue
        seen_genomes[cid].add(genome)
        d = row.get('depth_m', '')
        if d:
            depths[cid].append(float(d))

for c in clusters:
    ds = depths[c['full_id']]
    c['mean_depth'] = sum(ds) / len(ds) if ds else None
    c['n_with_depth'] = len(ds)

for c in clusters:
    color = CAT_PALETTE[clusters.index(c) % len(CAT_PALETTE)]
    c['color'] = color

# ---- overlap detection (connected components within OVERLAP_KM) ----
n = len(clusters)
adj = [[] for _ in range(n)]
for i in range(n):
    for j in range(i + 1, n):
        d = haversine_km(clusters[i]['lat'], clusters[i]['lon'], clusters[j]['lat'], clusters[j]['lon'])
        if d < OVERLAP_KM:
            adj[i].append(j)
            adj[j].append(i)

visited = [False] * n
groups = []
for i in range(n):
    if visited[i]:
        continue
    stack, comp = [i], []
    visited[i] = True
    while stack:
        node = stack.pop()
        comp.append(node)
        for nb in adj[node]:
            if not visited[nb]:
                visited[nb] = True
                stack.append(nb)
    groups.append(comp)

# assign display (label marker) position -- offset in a small ring for groups >1, unchanged for singletons
for comp in groups:
    if len(comp) == 1:
        i = comp[0]
        clusters[i]['disp_lat'] = clusters[i]['lat']
        clusters[i]['disp_lon'] = clusters[i]['lon']
        clusters[i]['needs_leader'] = False
    else:
        clat = sum(clusters[i]['lat'] for i in comp) / len(comp)
        clon = sum(clusters[i]['lon'] for i in comp) / len(comp)
        for k, i in enumerate(comp):
            angle = 2 * math.pi * k / len(comp) - math.pi / 2
            clusters[i]['disp_lat'] = clat + OFFSET_DEG * math.sin(angle)
            clusters[i]['disp_lon'] = clon + OFFSET_DEG * math.cos(angle) / math.cos(math.radians(clat))
            clusters[i]['needs_leader'] = True

# ---- figure ----
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})

fig = plt.figure(figsize=(18, 9), dpi=300)
MAP_LEFT = 0.29
ax = fig.add_axes([MAP_LEFT, 0.08, 0.99 - MAP_LEFT, 0.78], projection=ccrs.PlateCarree())
ax.set_global()
ax.add_feature(cfeature.OCEAN, facecolor='#E4EEEC', zorder=0)
ax.add_feature(cfeature.LAND, facecolor='#E9E4D6', edgecolor='#B9AF98', linewidth=0.5, zorder=1)
ax.gridlines(draw_labels=False, linewidth=0.4, color='#C7D0CB', linestyle=(0, (1, 3)), zorder=1)
for spine in ax.spines.values():
    spine.set_edgecolor('#3A4442')
    spine.set_linewidth(0.9)

for c in clusters:
    size = 60 + 55 * math.sqrt(c['n_genomes'] / max(cl['n_genomes'] for cl in clusters))
    if c['needs_leader']:
        ax.plot([c['lon'], c['disp_lon']], [c['lat'], c['disp_lat']], transform=ccrs.Geodetic(),
                 color=c['color'], linewidth=1.1, zorder=2, alpha=0.8)
        ax.scatter(c['lon'], c['lat'], transform=ccrs.PlateCarree(), s=26, color=c['color'],
                   edgecolor='white', linewidth=0.8, zorder=3)
        ax.scatter(c['disp_lon'], c['disp_lat'], transform=ccrs.PlateCarree(), s=size, color=c['color'],
                   edgecolor='#20302C', linewidth=1.0, zorder=4)
    else:
        ax.scatter(c['lon'], c['lat'], transform=ccrs.PlateCarree(), s=size, color=c['color'],
                   edgecolor='#20302C', linewidth=1.0, zorder=4)

    label_text = c['genus']
    ax.annotate(label_text, xy=(c['disp_lon'], c['disp_lat']), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                xytext=(6, 6), textcoords='offset points', fontsize=8.3, fontweight='bold',
                color='#20302C', style='italic', zorder=5,
                path_effects=[pe.withStroke(linewidth=2.2, foreground='white')])

fig.suptitle('PhaC Regional Specialists: 70%-Identity Clusters Concentrated in One Place',
             fontsize=17, fontweight='bold', x=(MAP_LEFT + 0.99) / 2, y=0.975)
fig.text((MAP_LEFT + 0.99) / 2, 0.905,
          'Selected from clusters with geo_mean_resultant_length > 0.85 (a genuinely concentrated, not cosmopolitan, footprint)\n'
          'and ≥8 genomes with resolved depth — position = circular-mean centroid of member genomes',
          ha='center', fontsize=9.8, color='#5B6E70')

# ---- legend: color, "genome" (dominant genus), average depth ----
legend_lines = []
for c in sorted(clusters, key=lambda c: (c['mean_depth'] is None, c['mean_depth'] or 0)):
    depth_str = f"{c['mean_depth']:.0f} m avg depth (n={c['n_with_depth']})" if c['mean_depth'] is not None else 'no depth data'
    legend_lines.append((c['color'], c['genus'], depth_str, c['n_genomes']))

legend_x, legend_y0, dy = 0.03, 0.10, 0.058
fig.text(legend_x, legend_y0 + dy * len(legend_lines) + 0.03, 'Genome\n(dominant genus)', fontsize=11.5, fontweight='bold', color='#20302C', va='bottom')
fig.text(legend_x + 0.145, legend_y0 + dy * len(legend_lines) + 0.03, 'Avg. depth', fontsize=11.5, fontweight='bold', color='#20302C', va='bottom')
for i, (color, genus, depth_str, n_genomes) in enumerate(legend_lines):
    y = legend_y0 + dy * (len(legend_lines) - 1 - i)
    fig.patches.append(plt.Rectangle((legend_x, y - 0.012), 0.018, 0.024, transform=fig.transFigure,
                                      facecolor=color, edgecolor='#20302C', linewidth=0.7, zorder=10))
    fig.text(legend_x + 0.028, y, f'{genus}', fontsize=11, style='italic', fontweight='bold', va='center')
    fig.text(legend_x + 0.145, y, f'{depth_str}', fontsize=10, color='#5B6E70', va='center')

footnote = (
    "Circle size ~ genome count (sqrt-scaled). Sydney, Australia pair (SHLQ01) sit ~1 km apart — shown as small true-position\n"
    "dots connected by a thin leader line to offset, larger label markers, rather than overlapping illegibly. Average depth is the\n"
    "mean depth_m across distinct genomes with resolved depth in that cluster (missing values dropped, not imputed as zero)."
)
fig.text((MAP_LEFT + 0.99) / 2, 0.015, footnote, ha='center', va='bottom', fontsize=8.4, color='#5B6E70')

out_path = OUT / 'phaC_regional_specialists.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phaC_regional_specialists.pdf', facecolor='white')
print('saved pdf too')

for c in clusters:
    print(f"{c['id']}  {c['genus']:<22}  lat={c['lat']:>7.2f} lon={c['lon']:>7.2f}  n={c['n_genomes']:>4}  "
          f"mean_depth={c['mean_depth']}  leader={c['needs_leader']}")
