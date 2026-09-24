"""Same radial dendrogram as plot_phac_taxonomy_radial.py (phaC 70%-
identity clusters found in exactly one host phylum -- dot size = cluster
count), with one addition: an outer ring of radial bars, one per phylum
leaf, showing the % of that phylum's phaC-positive GENOMES that carry
>=2 phaC copies (section 9's multi-copy analysis, not the single-phylum-
cluster count the dots already show -- a different question about the
same phyla, not a derived version of the first metric).

n_phaC per genome uses the same two corrections plot_phac_multicopy_genomes.py
applies (section 9, rebuilt 2026-09-24): decrement for genomes carrying
one of the 19 section-9.12-excluded no_hmm_triad_support target_ids, and
drop the 85 genomes the section-9.3 legitimacy audit called
likely_artifact outright. Recomputed fresh here rather than read back
from figures/phac_multicopy_by_phylum.tsv, since that file's own
pct_multicopy values can't be correctly re-aggregated into this script's
"N other phyla" grouping by averaging percentages (needs genome-level
weighting) -- easier and more correct to rebuild from the same raw join.

Usage:
    python figures/scripts/plot_phac_taxonomy_radial_multicopy.py

Output:
    figures/phac_taxonomy_radial_multicopy.png / .pdf
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _radial_dendrogram as rd
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import math

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
TOP_N_PHYLA = 18

# ---------------------------------------------------------------------
# 1. single-phylum cluster counts (exactly as plot_phac_taxonomy_radial.py)
# ---------------------------------------------------------------------
assignments = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        assignments[row[1]] = row[0]

cluster_phyla = defaultdict(set)
cluster_domain = {}
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if _phac_qc.is_bad(row['best_query']):
            continue
        cid = assignments.get(row['target_id'])
        if cid is None:
            continue
        if row.get('gtdb_phylum'):
            cluster_phyla[cid].add(row['gtdb_phylum'])
        if row.get('gtdb_domain') and cid not in cluster_domain:
            cluster_domain[cid] = row['gtdb_domain']

single_phylum_clusters = {cid: list(v)[0] for cid, v in cluster_phyla.items() if len(v) == 1}
n_single = len(single_phylum_clusters)
n_total = len(cluster_phyla)
print(f'{n_single}/{n_total} single-phylum clusters ({100*n_single/n_total:.1f}%)')

phylum_counts = defaultdict(int)
for cid, phylum in single_phylum_clusters.items():
    phylum_counts[phylum] += 1

ranked = sorted(phylum_counts.items(), key=lambda x: -x[1])
top_phyla = dict(ranked[:TOP_N_PHYLA])
other_phyla_names = {p for p, _ in ranked[TOP_N_PHYLA:]}
other_count = sum(n for p, n in ranked[TOP_N_PHYLA:])

paths = []
for phylum, n in top_phyla.items():
    dom = next((cluster_domain[cid] for cid, p in single_phylum_clusters.items() if p == phylum), 'Bacteria')
    paths.append(((dom, phylum), n))
other_label = f'{len(ranked) - TOP_N_PHYLA} other phyla'
if other_count:
    paths.append((('Bacteria', other_label), other_count))

root = rd.build_tree(paths)

# ---------------------------------------------------------------------
# 2. per-phylum %multicopy, current QC + section-9 corrections applied
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()

genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genomes[row['genome']] = {'n_phac': n, 'phylum': row['gtdb_phylum']}

genome_of_target = {}
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genome_of_target[row['target_id']] = row['genome']

bad_target_hits = defaultdict(int)
for tid, g in genome_of_target.items():
    if tid in _phac_qc.BAD_TARGET_IDS and g in genomes:
        bad_target_hits[g] += 1
for g, n_bad in bad_target_hits.items():
    genomes[g]['n_phac'] = max(0, genomes[g]['n_phac'] - n_bad)
genomes = {g: r for g, r in genomes.items() if r['n_phac'] > 0}

likely_artifact = set()
with open(OUT / 'phac_multicopy_legitimacy_audit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['legitimacy_call'] == 'likely_artifact':
            likely_artifact.add(row['genome'])
genomes = {g: r for g, r in genomes.items() if g not in likely_artifact}
print(f'{len(genomes):,} phaC-positive genomes after QC corrections (matches section 9)')


def pct_multicopy(phyla_names):
    recs = [r for r in genomes.values() if r['phylum'] in phyla_names]
    if not recs:
        return None, 0
    n_multi = sum(1 for r in recs if r['n_phac'] >= 2)
    return 100 * n_multi / len(recs), len(recs)


leaf_pct = {}
for phylum in top_phyla:
    pct, n = pct_multicopy({phylum})
    leaf_pct[phylum] = pct
    print(f'  {phylum:28s} pct_multicopy={pct:.1f}% (n={n:,} genomes)' if pct is not None else f'  {phylum}: no genomes')
if other_count:
    pct, n = pct_multicopy(other_phyla_names)
    leaf_pct[other_label] = pct
    print(f'  {other_label:28s} pct_multicopy={pct:.1f}% (n={n:,} genomes)' if pct is not None else f'  {other_label}: no genomes')

MAX_PCT_FOR_SCALE = 60.0  # observed range in section 9 tops out ~53%; round up for headroom

# ---------------------------------------------------------------------
# 3. render: base dendrogram, then a radial bar per leaf at its own angle
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(14, 14), dpi=300)

RADIUS_STEP = 1.3
RING_GAP = 0.18       # gap between leaf dot and start of its multicopy bar (units of RADIUS_STEP)
RING_MAX_LEN = 0.55   # max bar length at 100% (i.e. MAX_PCT_FOR_SCALE) (units of RADIUS_STEP)

angle_of, leaf_radius = rd.draw_radial_dendrogram(
    ax, root, radius_step=RADIUS_STEP, edge_color='#9AA5A0', edge_width=1.2,
    label_levels=(2,), label_fontsize=10.5,
    leaf_size_fn=lambda path, count: 30 + 260 * (count / max(top_phyla.values())),
    leaf_color_fn=lambda path, count: '#7A7A7A' if 'other' in path[-1] else '#C9622D',
    label_radius_offset=0.12 + RING_GAP + RING_MAX_LEN,  # push labels past the new ring
)

# draw_radial_dendrogram's own xlim/ylim margin (1.35x leaf radius) doesn't
# know about this script's extra ring + pushed-out labels, so the
# almost-straight-up Pseudomonadota label ran into the figure-level
# subtitle text (confirmed live) -- widen the axes' own limits to give
# the pushed-out labels the same margin the base dendrogram gives its own.
label_lim = leaf_radius + (0.12 + RING_GAP + RING_MAX_LEN + 0.55) * RADIUS_STEP
ax.set_xlim(-label_lim, label_lim)
ax.set_ylim(-label_lim, label_lim)


def leaf_nodes(node, path=()):
    count, children = node
    if not children:
        yield path, node
    for label, child in children.items():
        yield from leaf_nodes(child, path + (label,))


BAR_COLOR = '#1E6E7A'
for path, node in leaf_nodes(root):
    name = path[-1]
    pct = leaf_pct.get(name)
    if pct is None:
        continue
    angle = angle_of[id(node)]
    r0 = leaf_radius + RING_GAP * RADIUS_STEP
    r1 = r0 + RING_MAX_LEN * RADIUS_STEP * (pct / MAX_PCT_FOR_SCALE)
    x0, y0 = r0 * math.cos(angle), r0 * math.sin(angle)
    x1, y1 = r1 * math.cos(angle), r1 * math.sin(angle)
    ax.plot([x0, x1], [y0, y1], color=BAR_COLOR, linewidth=5.5, solid_capstyle='round', zorder=2, alpha=0.85)

# reference guide arcs at 20%/40%/60% multicopy, so bar lengths are readable.
# Label angle picked to land in clear whitespace, not on top of a branch --
# straight along angle=0 (used on the first pass) ran the labels right
# through the Pseudomonadota branch elbow.
REF_LABEL_ANGLE = -0.62  # radians, lower-right quadrant, clear of every leaf branch
for pct_ref in (20, 40, 60):
    r_ref = leaf_radius + RING_GAP * RADIUS_STEP + RING_MAX_LEN * RADIUS_STEP * (pct_ref / MAX_PCT_FOR_SCALE)
    n = 400
    xs = [r_ref * math.cos(2 * math.pi * i / n) for i in range(n + 1)]
    ys = [r_ref * math.sin(2 * math.pi * i / n) for i in range(n + 1)]
    ax.plot(xs, ys, color='#D8DDD9', linewidth=0.7, zorder=0, linestyle=':')
    lx, ly = r_ref * math.cos(REF_LABEL_ANGLE), r_ref * math.sin(REF_LABEL_ANGLE)
    ax.text(lx, ly, f'{pct_ref}%', fontsize=7.5, color='#8B958F', ha='center', va='bottom', zorder=0)

legend_handles = [
    mlines.Line2D([], [], marker='o', color='none', markerfacecolor='#C9622D', markersize=10,
                  label='Dot size: # host-specific (single-phylum) phaC clusters'),
    mlines.Line2D([], [], color=BAR_COLOR, linewidth=5.5, solid_capstyle='round',
                  label=f'Bar length: % of that phylum\'s phaC-positive genomes with >=2 copies (0-{MAX_PCT_FOR_SCALE:.0f}% scale)'),
]
ax.legend(handles=legend_handles, loc='upper left', fontsize=9.5, frameon=False, ncol=1)

fig.suptitle('phaC Host-Specific Clusters by Phylum, with Multi-Copy Prevalence', fontsize=17, fontweight='bold', y=0.97)
fig.text(0.5, 0.925,
          f'{n_single:,} of {n_total:,} phaC clusters ({100*n_single/n_total:.1f}%) are found in exactly one host phylum — top {TOP_N_PHYLA} shown individually. '
          'Outer teal bars: % of phaC-positive genomes in that phylum carrying >=2 copies (section 9).',
          ha='center', fontsize=10, color='#5B6E70', wrap=True)
fig.text(0.5, 0.015, 'Dot size = host-specific (single-phylum) phaC clusters. Bar length = genome-level multi-copy rate (an independent metric, same phyla). Innermost ring = domain.',
          ha='center', fontsize=8.7, color='#5B6E70')

out_path = OUT / 'phac_taxonomy_radial_multicopy.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_taxonomy_radial_multicopy.pdf', facecolor='white')
