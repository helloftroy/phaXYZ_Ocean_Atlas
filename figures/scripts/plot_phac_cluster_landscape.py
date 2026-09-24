"""Overview figure: the whole phaC_cluster0.7 (70%-identity) cluster
landscape, one point per cluster -- a bird's-eye view of the same
"species"-as-cluster framing section 7's rank-abundance analysis uses,
here spread across recurrence (how many genomes/how widely sampled) and
taxonomic breadth (how many different lineages independently carry it)
at once, rather than one dimension at a time.

x = log10(n_genomes carrying the cluster) -- abundance/recurrence.
y = n_distinct sampling locations carrying the cluster (asinh-scaled:
    log-like spread for the long tail, but doesn't collapse the dense 1-5
    range the way a log axis would) -- geographic recurrence. A location
    is a rounded (lat, lon) pair, not a sample_id or study_id: checked
    directly that sample_id over-counts (1,154 of 3,086 distinct
    locations in this dataset have >1 sample_id at the identical
    coordinate -- repeat visits, multiple depths/size-fractions at the
    same station) and study_id can under- or over-count in either
    direction (one study, e.g. TARA Oceans, spans many locations;
    conversely several independent studies can concentrate on one small
    region -- UBA10347 below is exactly that case, 14 studies but still
    geographically concentrated per its own geo_mean_resultant_length).
color = n_distinct GTDB classes among the cluster's genomes, bucketed
    (1/2/3/4+) rather than a continuous scale: nearly every cluster
    (97.6%) is confined to a single phylum, confirming the flat-color
    picture that level would give, but one level down (class) still
    shows real variation in the tail (95.8% single-class, but 148/20/8
    clusters span 2/3/4+ classes) -- informative without one level being
    so uniform it wastes the whole color channel on noise.

Cluster membership and genome counts are recomputed fresh from
phaC_cluster0.7_cluster.tsv + phaC_all_genomes_from_nr100_clusters.tsv
+ genome_family_matrix.tsv, filtered through the current
_phac_qc.load_bad_targets() (which now also excludes the 19 section-9.12
structurally-confirmed non-phaC targets) -- not read from the older
phaC_cluster0.7_cluster_ecology.tsv summary file, which predates that
exclusion and also only ever recorded each cluster's top-5 genera/phyla/
studies rather than true distinct counts.

Two groups highlighted with a dashed bounding box, per the two specific
cases this project expanded on earlier: the 16 regional-specialist
clusters (figures/phaC_regional_specialists.png, section 5 --
geographically concentrated, mostly single-location too) and the 3
Sulfitobacter clusters with >200 genomes each that were found to be
genuinely global instead (figures/phaC_sulfitobacter_overlay.png,
section 5.2) -- the same cluster IDs those two figures already used,
picked out here again by the identical selection logic rather than a
hardcoded list where that logic already exists in a script.

Usage:
    python figures/scripts/plot_phac_cluster_landscape.py

Output:
    figures/phac_cluster_landscape.png / .pdf
    figures/phac_cluster_landscape.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

REGIONAL_SPECIALIST_SUFFIXES = {
    '000004410215', '000199571761', '000234626934', '000033884939', '000088671462', '000036533161',
    '000077443878', '000226102035', '000027486053', '000159242036', '000128699032', '000183127334',
    '000236050034', '000159611932', '000192647371', '000113256229',
}
SULFITOBACTER_MIN_GENOMES = 200
# figures/scripts/plot_genus_global_overlay.py's >200-genome Sulfitobacter
# selection was run against the older phaC_cluster0.7_cluster_ecology.tsv;
# rederiving that same >200 threshold fresh here (current QC applied)
# drops the third cluster to 143 genomes (data attrition from QC
# corrections made since that figure), so the 3 cluster IDs are pinned
# explicitly instead -- these are the same 3 clusters section 5.2's
# "genuinely global" finding is about, not a fresh threshold re-derivation.
SULFITOBACTER_CLUSTER_IDS = {
    'OMDBv2.0_AA_G_NR100_000018755333', 'OMDBv2.0_AA_G_NR100_000045457530', 'OMDBv2.0_AA_G_NR100_000126950826',
}

# ---------------------------------------------------------------------
# build cluster -> {genomes} fresh, current QC applied
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()

genome_of = {}
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genome_of[row['target_id']] = row['genome']

gmeta = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        gmeta[row['genome']] = row

cluster_genomes = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        cid, tid = row[0], row[1]
        if tid in bad_targets:
            continue
        g = genome_of.get(tid)
        if g:
            cluster_genomes[cid].add(g)

print(f'{len(cluster_genomes):,} clusters with >=1 QC-passing genome')

def location_of(r):
    """Rounded (lat, lon) as a proxy for "one sampling location" -- distinct
    from sample_id, which over-counts: 1,154 of 3,086 distinct locations in
    this dataset have >1 sample_id at the identical coordinate (repeat
    visits, multiple depths/size-fractions at the same station), confirmed
    directly against genome_family_matrix.tsv before picking this over the
    simpler sample_id count."""
    if not r['latitude_degN'] or not r['longitude_degE']:
        return None
    return (round(float(r['latitude_degN']), 3), round(float(r['longitude_degE']), 3))


clusters = []
for cid, genomes in cluster_genomes.items():
    known = [gmeta[g] for g in genomes if g in gmeta]
    if not known:
        continue
    n_genomes = len(genomes)
    n_locations = len({loc for r in known if (loc := location_of(r)) is not None})
    n_classes = len(set(r['gtdb_class'] for r in known))
    top_genus = max(set(r['gtdb_genus'] for r in known), key=lambda g: sum(1 for r in known if r['gtdb_genus'] == g))
    clusters.append({'cluster_id': cid, 'n_genomes': n_genomes, 'n_locations': n_locations,
                      'n_classes': n_classes, 'top_genus': top_genus})

clusters.sort(key=lambda c: -c['n_genomes'])
out_tsv = OUT / 'phac_cluster_landscape.tsv'
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(clusters[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(clusters)
print(f'wrote {out_tsv} ({len(clusters):,} clusters)')

regional_ids = {c['cluster_id'] for c in clusters if c['cluster_id'][-12:] in REGIONAL_SPECIALIST_SUFFIXES}
sulfitobacter_ids = {c['cluster_id'] for c in clusters if c['cluster_id'] in SULFITOBACTER_CLUSTER_IDS}
print(f'{len(regional_ids)}/16 regional-specialist clusters found, {len(sulfitobacter_ids)}/3 Sulfitobacter global clusters found')
for cid in sulfitobacter_ids:
    c = next(c for c in clusters if c['cluster_id'] == cid)
    if c['n_genomes'] <= SULFITOBACTER_MIN_GENOMES:
        print(f'  note: {cid} now has only {c["n_genomes"]} genomes post-QC (below the original >200 threshold)')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
CLASS_COLOR = {1: '#B7BDB8', 2: '#4A7FB5', 3: '#D98E04', 4: '#9E3B3B'}
CLASS_LABEL = {1: '1 class', 2: '2 classes', 3: '3 classes', 4: '4+ classes'}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5})
fig, ax = plt.subplots(figsize=(12, 9), dpi=300)

def bucket_of(n_classes):
    return min(n_classes, 4)


# jitter: n_genomes and n_locations are both small integers, so unjittered
# points overlap heavily into a sparse-looking grid (confirmed live on
# the first render -- most of the 4,340 clusters were invisible, stacked
# under a handful of grid points) -- a small multiplicative/additive
# jitter spreads them into a readable cloud without moving any point
# more than a fraction of the gap to its neighboring integer.
rng = np.random.default_rng(0)
for bucket in (1, 2, 3, 4):
    pts = [c for c in clusters if bucket_of(c['n_classes']) == bucket]
    xs = np.array([c['n_genomes'] for c in pts], dtype=float) * rng.uniform(0.85, 1.18, len(pts))
    ys = np.array([c['n_locations'] for c in pts], dtype=float) + rng.uniform(-0.32, 0.32, len(pts))
    zorder = 2 + bucket
    ax.scatter(xs, ys, s=26 if bucket > 1 else 14, color=CLASS_COLOR[bucket],
               alpha=0.9 if bucket > 1 else 0.3, linewidth=0, zorder=zorder, label=CLASS_LABEL[bucket])

ax.set_xscale('log')
ax.set_yscale('asinh')
# asinh is symmetric about 0 by default, which -- since every real value
# here is >=1 -- pulled in a whole unused negative-axis half full of
# ticks (-10^0, -10^-1, ...) that confirmed live looked like a rendering
# bug on the first pass. Explicit ylim + explicit positive-only ticks
# (plain integers, not scientific notation) fixes it.
ax.set_ylim(0.55, 350)
YTICKS = [1, 2, 3, 5, 10, 20, 50, 100, 200, 300]
ax.set_yticks(YTICKS)
ax.set_yticklabels([str(t) for t in YTICKS])
ax.set_xlabel('Genomes carrying the cluster (log scale)')
ax.set_ylabel('Distinct sampling locations carrying the cluster (asinh scale)')
ax.set_title('The phaC cluster landscape: recurrence vs. geographic spread vs. taxonomic breadth',
              fontsize=15, fontweight='bold')
ax.grid(True, which='major', color='#E4E8E5', linewidth=0.6, zorder=0)
ax.grid(True, which='minor', color='#F0F2F0', linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


def draw_box(ids, color, label, label_xy, pad_x=1.35, pad_y_add=0.6):
    pts = [c for c in clusters if c['cluster_id'] in ids]
    x0, x1 = min(c['n_genomes'] for c in pts) / pad_x, max(c['n_genomes'] for c in pts) * pad_x
    y0, y1 = max(0.5, min(c['n_locations'] for c in pts) - pad_y_add), max(c['n_locations'] for c in pts) + pad_y_add
    rect = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle='round,pad=0,rounding_size=0',
                            linewidth=1.8, edgecolor=color, facecolor='none', linestyle='--', zorder=10)
    ax.add_patch(rect)
    ax.annotate(label, xy=(x1, y1), xytext=label_xy, textcoords='data', fontsize=10, color=color, fontweight='bold',
                ha='left', arrowprops=dict(arrowstyle='-', color=color, lw=1.2, shrinkA=0, shrinkB=4))
    # ring the exact points, not just the box -- the box necessarily also
    # encloses unrelated clusters that happen to share similar coordinates
    ax.scatter([c['n_genomes'] for c in pts], [c['n_locations'] for c in pts], s=90, facecolor='none',
               edgecolor=color, linewidth=1.4, zorder=11)
    return pts


draw_box(regional_ids, '#7A5FA0', '16 regional specialists\n(single region, but sometimes many nearby\nstations -- section 5)', (2.2, 130))
draw_box(sulfitobacter_ids, '#B33951', 'the 3 truly global\nSulfitobacter clusters (section 5.2)', (9, 220))

class_handles = [mpatches.Patch(color=CLASS_COLOR[b], label=CLASS_LABEL[b]) for b in (1, 2, 3, 4)]
leg1 = ax.legend(handles=class_handles, loc='lower right', title='Color: # distinct GTDB classes', fontsize=9.5,
                   title_fontsize=10, frameon=True, facecolor='white', edgecolor='#D8DDD9')
ax.add_artist(leg1)

fig.text(0.5, 0.005, f'{len(clusters):,} phaC_cluster0.7 clusters (70%-identity groups), one point each, current QC applied '
                       f'(figures/scripts/_phac_qc.py, including section 9.12\'s no_hmm_triad_support exclusion).',
          ha='center', fontsize=8.7, color='#5B6E70')
fig.tight_layout(rect=[0, 0.02, 1, 1])

out_path = OUT / 'phac_cluster_landscape.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_cluster_landscape.pdf', facecolor='white')
