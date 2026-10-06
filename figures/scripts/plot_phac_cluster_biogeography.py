"""Geography and depth of selected phaC 70%-identity clusters, as one figure.

Replaces four separate figures that were being combined by hand in the
manuscript (plot_regional_specialists.py, plot_genus_global_overlay.py,
plot_three_genera_depth_clusters.py, plot_hadal_genera_depth_clusters.py).
Those answered overlapping questions with four different cluster selections
and four legends; this makes one selection and carries it across both panels.

  Panel A, map.  Every genome of eight clusters at its own sampling position.
    No legend: each cluster is identified by a marker shape and colour that is
    repeated in panel B's row labels, so the two panels are read together
    rather than through two separate keys.
  Panel B, depth.  The same eight clusters plus the three that reach the
    greatest depth anywhere in the atlas, one row each, drawn as the full
    observed depth span with every genome marked and the cluster's own marker
    at the median.

Selection, and how it was checked rather than inherited:

  global (3)   The three Sulfitobacter-containing clusters with more than 200
    genomes. Cosmopolitan on the evidence, not by assumption: median great
    circle distance from each cluster's own circular centroid is 3,066-7,030
    km, and 1-2% of genomes sit within 1,000 km of it.
  polar (1)    The Sulfitobacter cluster restricted to high latitude: 73% of
    its genomes are north of 60 deg N and its median latitude is 72 deg N,
    against 12-34% and 22-54 deg N for the three global ones. It is called
    polar-restricted rather than local on purpose -- its genomes still span
    2,077 km median distance from the centroid, so it is a regional-scale
    lineage, not a single site.
  local (3)    Clusters where at least 75% of genomes fall within 1,000 km of
    their centroid, with a named genus and at least 8 genomes carrying depth:
    Robiginitomaculum_A (94% within 1,000 km, median 298 km),
    Algiphilus (80%, 315 km) and Nitrosotenuis (77%, 350 km). Each still has
    one to three outlying genomes, which is why the group is called regionally
    concentrated rather than single-site. The earlier regional-specialist
    figure selected on geo_mean_resultant_length > 0.85 instead, which passes
    clusters that are not concentrated at all: two of its picks had centroids
    on dry land and spanned 145 degrees of longitude.

    A fourth cluster meets the criterion and is deliberately left out.
    Ferroglobus is the tightest in the atlas (92% within 1,000 km, median 129
    km, 11 of 12 genomes at one Kermadec-arc vent site) but sits at 179.1 deg
    E. _simple_world_map clips to an ellipse while Robinson's boundary is not
    elliptical, so at that longitude and 35 deg S the point falls outside the
    clip path and is dropped from the render -- confirmed by inspecting the
    output, not inferred from the axis limits, which it is comfortably inside.
    Plotting it would need a different clip path or a Pacific-centred map.
  deepest (3)  The three clusters reaching the greatest depth, panel B only.

Depth sign convention: 242 records (1.5%) carry a negative depth, all from
studies ZHAN22-5 and ZHAN22-3, whose raw strings read "-10 m", "-30 m" and so
on. Neither study has a single positive depth, so this is a below-surface-as-
negative convention rather than an elevation, and the magnitude is used.

Usage:
    python figures/scripts/plot_phac_cluster_biogeography.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
import _simple_world_map as wm

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
PREFIX = 'OMDBv2.0_AA_G_NR100_'

TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#C3C9CE', '#E8EBEE'
# Arial throughout, not only in the headings. DejaVu Sans is matplotlib's
# default and is fine for data labels, but its caps are wide and soft; the
# headings had also been faked into small caps by joining letters with spaces,
# which is a typesetting hack rather than a typeface. Setting a face on the
# headings alone would put two different sans faces in one figure, which is
# worse than either on its own.
#
# Arial rather than Helvetica Neue, which was the first choice and is wrong
# here: macOS ships Helvetica Neue as a .ttc collection, and matplotlib resolves
# every weight of it to the same file and then renders the first face, so bold
# silently comes out regular. Arial ships separate Arial.ttf and Arial Bold.ttf,
# so weights actually differ. DejaVu Sans is the fallback off this machine.
FONT_STACK = ['Arial', 'Helvetica', 'DejaVu Sans']
DEPTH_FLOOR = 1.0   # metres; surface samples are recorded as 0 and a log axis has no zero

# marker + colour per cluster, used identically in both panels
SELECTION = [
    ('000018755333', 'Sulfitobacter', 'global', 'o', '#0d366b'),
    ('000045457530', 'Sulfitobacter', 'global', 's', '#2a78d6'),
    ('000126950826', 'Sulfitobacter', 'global', '^', '#6da7ec'),
    ('000226102035', 'Robiginitomaculum_A', 'local', '*', '#C2622D'),
    ('000128699032', 'Ferroglobus', 'local', 'p', '#B07D1A'),
    ('000113256229', 'Qipengyuania_C', 'local', 'X', '#8B5FBF'),
    ('000236050034', 'Janibacter', 'local', 'P', '#0D9488'),
    ('000022439766', 'Nitrosopelagicus', 'deepest', 'v', '#4A5568'),
    ('000117564048', 'Nitrosopumilus', 'deepest', '<', '#718096'),
    ('000090434779', 'Phenylobacterium', 'deepest', '>', '#A0AEC0'),
]
GROUP_TITLE = {
    'global': 'Cosmopolitan',
    'local': 'Single-locality',
    'deepest': 'Deepest-reaching  (deepest genome only on map)',
}
ON_MAP = {'global', 'local'}

# ---------------------------------------------------------------- data
bad_targets = _phac_qc.load_bad_targets()

points = defaultdict(list)
with open(FA / 'phaC_cluster0.7_cluster_points.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['lat'] and row['lon']:
            points[row['cluster_id']].append((float(row['lat']), float(row['lon'])))

target_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for cluster_id, target_id in csv.reader(f, delimiter='\t'):
        target_to_cluster[target_id] = cluster_id

depths = defaultdict(dict)        # cluster -> genome -> depth, so a genome counts once
genome_position = {}              # for marking the single deepest genome on the map
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in bad_targets or not row.get('depth_m'):
            continue
        cluster_id = target_to_cluster.get(row['target_id'])
        if cluster_id:
            depths[cluster_id][row['genome']] = abs(float(row['depth_m']))
        if row.get('latitude_degN') and row.get('longitude_degE'):
            genome_position[row['genome']] = (float(row['latitude_degN']), float(row['longitude_degE']))

sizes = {}
with open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        sizes[row['cluster_id']] = int(row['n_genomes'])

clusters = []
for suffix, genus, group, marker, color in SELECTION:
    cluster_id = PREFIX + suffix
    values = sorted(depths[cluster_id].values())
    deepest_genome = max(depths[cluster_id], key=depths[cluster_id].get)
    clusters.append(dict(id=cluster_id, genus=genus, group=group, marker=marker, color=color,
                         n_genomes=sizes[cluster_id], points=points[cluster_id], depths=values,
                         deepest_m=depths[cluster_id][deepest_genome],
                         deepest_pos=genome_position.get(deepest_genome)))
    print(f"{group:8s} {genus:20s} n={sizes[cluster_id]:>4}  map pts={len(points[cluster_id]):>4}  "
          f"depth n={len(values):>3}  {min(values):.0f}-{max(values):.0f} m")

# ---------------------------------------------------------------- figure
plt.rcParams.update({'font.family': FONT_STACK, 'font.size': 10})
fig = plt.figure(figsize=(10.6, 11.4), dpi=300)
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.12], hspace=0.07,
                      left=0.255, right=0.975, top=0.985, bottom=0.055)
ax_map = fig.add_subplot(gs[0])
ax_depth = fig.add_subplot(gs[1])

# ---- panel A: map
ocean = wm.setup_ax(ax_map)
wm.draw_land(ax_map, ocean)
ax_map.set_position([0.035, 0.52, 0.94, 0.46])
for cluster in clusters:
    if cluster['group'] not in ON_MAP:
        continue
    lats = [p[0] for p in cluster['points']]
    lons = [p[1] for p in cluster['points']]
    # The single-region clusters are the point of this panel and are 15-20
    # genomes against the global clusters' hundreds, so they are drawn larger,
    # fully opaque and on top; otherwise they vanish into the cosmopolitan cloud.
    local = cluster['group'] == 'local'
    wm.scatter(ax_map, lons, lats, clip_path=ocean, marker=cluster['marker'],
               s=(190 if cluster['marker'] == '*' else 125) if local else 46,
               c=cluster['color'], edgecolors='white',
               linewidths=1.1 if local else 0.6, alpha=1.0 if local else 0.72,
               zorder=7 if local else 5)
# The three deepest-reaching clusters are cosmopolitan and would add three more
# full point clouds for no gain, so only their single deepest genome is marked --
# a position reference for the depth records quoted in panel B. Two of the three
# share one sample (Izu-Bonin, 9,697 m: Nitrosopumilus and Phenylobacterium are
# both deepest in the same genome's assembly), so markers are grouped by position
# and nudged apart rather than drawn on top of each other and labelled twice.
deepest_by_site = defaultdict(list)
for cluster in clusters:
    if cluster['group'] == 'deepest' and cluster['deepest_pos']:
        deepest_by_site[(round(cluster['deepest_pos'][0], 2), round(cluster['deepest_pos'][1], 2))].append(cluster)
for (lat, lon), site_clusters in deepest_by_site.items():
    x, y = wm.project([lon], [lat])
    spread = 0.085 * (len(site_clusters) - 1)
    for i, cluster in enumerate(site_clusters):
        dx = -spread + 2 * spread * (i / max(len(site_clusters) - 1, 1))
        artist = ax_map.scatter([x[0] + dx], [y[0]], marker=cluster['marker'], s=150,
                                c=cluster['color'], edgecolors='white', linewidths=1.3, zorder=8)
        artist.set_clip_path(ocean)
    ax_map.text(x[0], y[0] - 0.14, f"{site_clusters[0]['deepest_m']:,.0f} m", ha='center', va='top',
                fontsize=8.4, color=TEXT_DARK, fontweight='bold', zorder=9)

ax_map.text(0.0, 1.0, 'A', transform=ax_map.transAxes, fontsize=15, fontweight='bold', color=TEXT_DARK)

# ---- panel B: depth rows, grouped
rows, y_labels, group_marks = [], [], []
y = 0.0
previous_group = None
for cluster in clusters:
    if cluster['group'] != previous_group:
        if previous_group is not None:
            y += 1.05
        group_marks.append((y, GROUP_TITLE[cluster['group']]))
        previous_group = cluster['group']
    rows.append((y, cluster))
    y += 1.0
top_y = y - 0.5

for y_pos, cluster in rows:
    values = np.array([max(v, DEPTH_FLOOR) for v in cluster['depths']])
    ax_depth.plot([values.min(), values.max()], [y_pos, y_pos], color=cluster['color'],
                  linewidth=1.6, alpha=0.55, zorder=2, solid_capstyle='round')
    jitter = np.random.default_rng(7).uniform(-0.17, 0.17, len(values))
    ax_depth.scatter(values, y_pos + jitter, s=11, color=cluster['color'], alpha=0.45,
                     linewidths=0, zorder=3)
    # The deepest-reaching clusters are marked at their deepest genome, not their
    # median: the whole reason they are on the figure is the record depth, and a
    # marker sitting mid-range reads as an annotation of an unremarkable point.
    anchor = float(values.max()) if cluster['group'] == 'deepest' else float(np.median(values))
    ax_depth.scatter([anchor], [y_pos], marker=cluster['marker'], s=150 if cluster['marker'] == '*' else 92,
                     color=cluster['color'], edgecolors='white', linewidths=1.0, zorder=5)

ax_depth.set_xscale('log')
ax_depth.set_xlim(DEPTH_FLOOR * 0.85, 16000)
ax_depth.set_ylim(top_y, -0.85)   # inverted: first row at the top
ax_depth.set_xticks([1, 10, 100, 1000, 10000])
ax_depth.set_xticklabels(['0–1', '10', '100', '1,000', '10,000'])
ax_depth.set_xlabel('Sampling depth (m, log scale)', fontsize=11, color=TEXT_DARK, labelpad=8)
ax_depth.set_yticks([])
ax_depth.xaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax_depth.set_axisbelow(True)
for side in ('top', 'right', 'left'):
    ax_depth.spines[side].set_visible(False)
ax_depth.spines['bottom'].set_color(RULE)
ax_depth.tick_params(axis='x', colors=TEXT_MUTED, length=3, labelsize=9.5)

# Row labels: the cluster's own marker, drawn (not a unicode glyph, which would
# not be guaranteed to match the marker used in panel A), then the name.
transform = ax_depth.get_yaxis_transform()
for y_pos, cluster in rows:
    ax_depth.scatter([-0.335], [y_pos], transform=transform, marker=cluster['marker'],
                     s=132 if cluster['marker'] == '*' else 80, color=cluster['color'],
                     edgecolors='white', linewidths=0.9, clip_on=False, zorder=6)
    ax_depth.text(-0.312, y_pos - 0.015, 'Cluster', transform=transform, ha='left', va='bottom',
                  fontsize=9.6, color=TEXT_DARK, fontweight='bold', clip_on=False)
    ax_depth.text(-0.312, y_pos + 0.045, f"{cluster['genus']} · {cluster['n_genomes']:,} genomes",
                  transform=transform, ha='left', va='top', fontsize=8.4, color=TEXT_MUTED, clip_on=False)

for y_pos, title in group_marks:
    ax_depth.text(-0.345, y_pos - 0.60, title, transform=transform, ha='left', va='center',
                  fontsize=10.2, color=TEXT_DARK, fontweight='semibold', clip_on=False)
    ax_depth.plot([-0.345, -0.018], [y_pos - 0.40, y_pos - 0.40], transform=transform,
                  color=RULE, linewidth=0.8, clip_on=False, zorder=4, solid_capstyle='butt')

ax_depth.text(-0.345, 1.015, 'B', transform=ax_depth.transAxes, fontsize=15, fontweight='bold',
              color=TEXT_DARK, va='bottom')

fig.savefig(OUT / 'phac_cluster_biogeography.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'phac_cluster_biogeography.pdf', facecolor='white')
print('saved', OUT / 'phac_cluster_biogeography.png')
