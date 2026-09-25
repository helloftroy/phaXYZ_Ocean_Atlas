"""Same Phylum->Class->Order dendrogram as plot_phac_taxonomy_radial_detailed.py
(all 38 phyla, top 3 classes/phylum, top 2 orders/class, 104 leaves,
outer colored band = phylum), with an outer ring of radial bars added
past the phylum band: % of phaC-positive genomes matching that EXACT
leaf's taxonomy (order, or the folded "N smaller orders"/"(mixed)"
bucket it represents) carrying >=2 phaC copies. 104 bars instead of the
18 the phylum-only version (plot_phac_taxonomy_radial_multicopy.py) has
-- the whole reason for building this against the detailed tree instead:
order-level resolution shows real within-phylum variation the phylum-
level version necessarily averages away.

n_phaC per genome uses the same two corrections section 9's multicopy
analysis applies (decrement for genomes carrying one of the 19 section-
9.12-excluded no_hmm_triad_support target_ids; drop the 85 likely_artifact
genomes outright). Genome-level (gtdb_phylum, gtdb_class, gtdb_order) is
looked up directly per genome -- NOT via the cluster-level class/order
assignment the tree itself is built from (which is single-valued per
cluster and only used for the tree's own shape) -- so a leaf's bar
reflects every genome actually in that taxon, not just genomes that
happen to carry a cluster assigned to it.

Usage:
    python figures/scripts/plot_phac_taxonomy_radial_detailed_multicopy.py

Output:
    figures/phac_taxonomy_radial_detailed_multicopy.png / .pdf
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _radial_dendrogram as rd
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.lines as mlines
from matplotlib.patches import Wedge
import colorsys

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
TOP_K_CLASSES = 3
TOP_K_ORDERS = 2
OTHER_COLOR = '#7A7A7A'

# ---------------------------------------------------------------------
# 1. tree structure -- identical logic to plot_phac_taxonomy_radial_detailed.py,
#    but also records, per leaf, exactly which (class, order) combinations
#    it represents, so genome-level stats can be matched to the same leaf.
# ---------------------------------------------------------------------
assignments = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        assignments[row[1]] = row[0]

cluster_phyla = defaultdict(set)
cluster_class = {}
cluster_order = {}
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
        if row.get('gtdb_class') and cid not in cluster_class:
            cluster_class[cid] = row['gtdb_class']
        if row.get('gtdb_order') and cid not in cluster_order:
            cluster_order[cid] = row['gtdb_order']
        if row.get('gtdb_domain') and cid not in cluster_domain:
            cluster_domain[cid] = row['gtdb_domain']

single_phylum = {cid: list(v)[0] for cid, v in cluster_phyla.items() if len(v) == 1}
n_single, n_total = len(single_phylum), len(cluster_phyla)
print(f'{n_single}/{n_total} single-phylum clusters ({100*n_single/n_total:.1f}%)')

phylum_totals = defaultdict(int)
phylum_domain = {}
class_counts_by_phylum = defaultdict(lambda: defaultdict(int))
order_counts_by_class = defaultdict(lambda: defaultdict(int))

for cid, phylum in single_phylum.items():
    phylum_totals[phylum] += 1
    phylum_domain[phylum] = cluster_domain.get(cid, 'Bacteria')
    cls = cluster_class.get(cid, 'Unknown class')
    class_counts_by_phylum[phylum][cls] += 1
    order = cluster_order.get(cid, 'Unknown order')
    order_counts_by_class[(phylum, cls)][order] += 1

n_phyla = len(phylum_totals)
print(f'{n_phyla} distinct phyla -- ALL shown individually (no "other phyla" bucket)')

paths = []
leaf_criteria = {}  # path -> ('order', phylum, cls, order) | ('other_orders', phylum, cls, {orders}) | ('other_classes', phylum, {classes})
for phylum, total in phylum_totals.items():
    dom = phylum_domain[phylum]
    classes_ranked = sorted(class_counts_by_phylum[phylum].items(), key=lambda x: -x[1])
    shown_classes = classes_ranked[:TOP_K_CLASSES]
    other_classes = classes_ranked[TOP_K_CLASSES:]

    for cls, cls_n in shown_classes:
        orders_ranked = sorted(order_counts_by_class[(phylum, cls)].items(), key=lambda x: -x[1])
        shown_orders = orders_ranked[:TOP_K_ORDERS]
        other_orders = orders_ranked[TOP_K_ORDERS:]
        for order, order_n in shown_orders:
            path = (dom, phylum, cls, order)
            paths.append((path, order_n))
            leaf_criteria[path] = ('order', phylum, cls, order)
        if other_orders:
            n_other = len(other_orders)
            other_n = sum(n for _, n in other_orders)
            path = (dom, phylum, cls, f'{n_other} smaller orders')
            paths.append((path, other_n))
            leaf_criteria[path] = ('other_orders', phylum, cls, {o for o, _ in other_orders})

    if other_classes:
        n_other_cls = len(other_classes)
        other_cls_n = sum(n for _, n in other_classes)
        path = (dom, phylum, f'{n_other_cls} smaller classes', '(mixed)')
        paths.append((path, other_cls_n))
        leaf_criteria[path] = ('other_classes', phylum, {c for c, _ in other_classes})

root = rd.build_tree(paths)
n_leaves = len(paths)
print(f'{n_leaves} leaves shown (order-level, or "smaller classes/orders" buckets)')

# ---------------------------------------------------------------------
# 2. per-leaf %multicopy, genome-level, current QC + section-9 corrections
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()

genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genomes[row['genome']] = {'n_phac': n, 'phylum': row['gtdb_phylum'], 'class': row['gtdb_class'], 'order': row['gtdb_order']}

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

# index genomes by phylum once for speed (104 leaves x scanning 31k genomes
# each would be slow; grouping by phylum first cuts the per-leaf scan a lot)
genomes_by_phylum = defaultdict(list)
for g, r in genomes.items():
    genomes_by_phylum[r['phylum']].append(r)


def matching_genomes(criterion):
    kind = criterion[0]
    phylum = criterion[1]
    pool = genomes_by_phylum.get(phylum, [])
    if kind == 'order':
        _, _, cls, order = criterion
        return [r for r in pool if r['class'] == cls and r['order'] == order]
    if kind == 'other_orders':
        _, _, cls, orders = criterion
        return [r for r in pool if r['class'] == cls and r['order'] in orders]
    if kind == 'other_classes':
        _, _, classes = criterion
        return [r for r in pool if r['class'] in classes]
    raise ValueError(criterion)


leaf_pct = {}
n_no_genomes = 0
for path, criterion in leaf_criteria.items():
    recs = matching_genomes(criterion)
    if not recs:
        n_no_genomes += 1
        continue
    n_multi = sum(1 for r in recs if r['n_phac'] >= 2)
    leaf_pct[path] = (100 * n_multi / len(recs), len(recs))
print(f'{len(leaf_pct)}/{n_leaves} leaves have >=1 matching phaC-positive genome ({n_no_genomes} have none post-QC)')

MAX_PCT_FOR_SCALE = 70.0

# ---------------------------------------------------------------------
# 3. render: base dendrogram + phylum ring (unchanged), then multicopy ring
# ---------------------------------------------------------------------
def gen_phylum_colors(phyla_sorted_by_size):
    golden_conjugate = 0.6180339887498949
    colors = {}
    hue = 0.15
    for i, phylum in enumerate(phyla_sorted_by_size):
        hue = (hue + golden_conjugate) % 1.0
        sat = 0.58 if i % 2 == 0 else 0.72
        val = 0.72 if i % 3 != 0 else 0.62
        r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
        colors[phylum] = (r, g, b)
    return colors


phyla_ranked = sorted(phylum_totals.items(), key=lambda x: -x[1])
PHYLUM_COLORS = gen_phylum_colors([p for p, _ in phyla_ranked])

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(22, 22), dpi=300)


def branch_color(path, count):
    if 'smaller' in path[-1] or 'smaller' in path[-2]:
        return OTHER_COLOR
    return PHYLUM_COLORS.get(path[1], OTHER_COLOR)


RADIUS_STEP = 1.05
MAX_DEPTH = 4  # domain/phylum/class/order -- deterministic, same fact the base detailed script relies on
MAX_R = MAX_DEPTH * RADIUS_STEP  # == what draw_radial_dendrogram will return as max_r; needed up front
RING_GAP = 0.05          # gap (units of MAX_R) between phylum band outer edge and multicopy ring start
RING_MAX_LEN = 0.20      # multicopy ring max length (units of MAX_R) at 100%-of-scale -- shortened
                          # from an original 0.34 per direct feedback that the ring read as too long/thin
LABEL_MARGIN = 0.06      # extra clearance (units of MAX_R) past the longest possible bar, before labels start
# label_radius_offset is in units of radius_step (draw_radial_dendrogram computes
# label_r = leaf_r + radius_step*offset -- NOT max_r*offset), so RING_GAP/RING_MAX_LEN/
# LABEL_MARGIN (defined in units of MAX_R above, to match how the ring itself is drawn)
# must be converted before being passed in -- a units mismatch here (max_r-relative
# fractions passed in raw as if they were radius_step-relative) was the actual bug
# behind the first render's labels landing inside the new ring instead of past it.
LABEL_RADIUS_OFFSET = (0.65 * RADIUS_STEP + (RING_GAP + RING_MAX_LEN + LABEL_MARGIN) * MAX_R) / RADIUS_STEP

angle_of, max_r = rd.draw_radial_dendrogram(
    ax, root, radius_step=RADIUS_STEP, edge_color='#D8DED7', edge_width=0.7,
    label_levels=(4,), label_fontsize=6.3, label_radius_offset=LABEL_RADIUS_OFFSET,
    leaf_size_fn=lambda path, count: 0,  # replaced by the thick tick marks below -- a low-count
                                          # leaf's dot (min size 6px) was confirmed live to be
                                          # nearly invisible among 104 leaves at this figure size
    leaf_color_fn=branch_color, edge_color_fn=branch_color,
    show_internal_dots=True, internal_dot_size=5, internal_dot_color='#DDE2DB',
)

def leaf_nodes(node, path=()):
    count, children = node
    if not children:
        yield path, node
    for label, child in children.items():
        yield from leaf_nodes(child, path + (label,))


# ---- leaf cluster-count marker: a short thick radial tick, not a dot --
# linewidth (not marker size) carries the count, so even a count near the
# bottom of the range stays visibly a solid mark rather than shrinking to
# a near-invisible point.
TICK_HALF_LEN = 0.045 * MAX_R
max_leaf_count = max(n for _, n in paths)
for path, node in leaf_nodes(root):
    count = node[0]
    angle = angle_of[id(node)]
    lw = 1.3 + 7.5 * (count / max_leaf_count)
    r0, r1 = max_r - TICK_HALF_LEN, max_r + TICK_HALF_LEN
    x0, y0 = r0 * math.cos(angle), r0 * math.sin(angle)
    x1, y1 = r1 * math.cos(angle), r1 * math.sin(angle)
    ax.plot([x0, x1], [y0, y1], color=branch_color(path, count), linewidth=lw, solid_capstyle='round', zorder=3)

# ---- outer colored clade ring, one wedge per phylum (depth 2) ----
ring_r0 = max_r * 1.05
ring_r1 = max_r * 1.11
phylum_nodes = []
for dom_label, dom_child in rd._sorted_children(root[1]):
    for label, child in rd._sorted_children(dom_child[1]):
        phylum_nodes.append((label, child))

for label, child in phylum_nodes:
    lo, hi = rd.leaf_angle_range(child, angle_of)
    pad = 0.003
    color = PHYLUM_COLORS.get(label, OTHER_COLOR)
    wedge = Wedge((0, 0), ring_r1, math.degrees(lo - pad), math.degrees(hi + pad),
                  width=ring_r1 - ring_r0, facecolor=color, alpha=0.85, edgecolor='white', linewidth=0.6, zorder=0)
    ax.add_patch(wedge)
    span_deg = math.degrees(hi - lo)
    if span_deg > 1.3:
        mid = (lo + hi) / 2
        lx, ly = rd._polar(mid, ring_r1 + max_r * 0.035)
        rot = math.degrees(mid)
        ha = 'left'
        if math.cos(mid) < 0:
            rot += 180
            ha = 'right'
        short_label = label if len(label) < 20 else label[:18] + '…'
        fontsize = 9.5 if span_deg > 6 else 6.8
        ax.text(lx, ly, short_label, rotation=rot, rotation_mode='anchor', ha=ha, va='center',
                fontsize=fontsize, fontweight='bold', color='#20302C', zorder=5,
                path_effects=[pe.withStroke(linewidth=2.2, foreground='white')])

# ---- outer multicopy ring: one radial bar per leaf, starting past the phylum band ----
BAR_COLOR = '#1E6E7A'
mc_r0 = ring_r1 + RING_GAP * max_r
for path, node in leaf_nodes(root):
    pct_n = leaf_pct.get(path)
    if pct_n is None:
        continue
    pct, n = pct_n
    angle = angle_of[id(node)]
    r1 = mc_r0 + RING_MAX_LEN * max_r * (pct / MAX_PCT_FOR_SCALE)
    x0, y0 = mc_r0 * math.cos(angle), mc_r0 * math.sin(angle)
    x1, y1 = r1 * math.cos(angle), r1 * math.sin(angle)
    ax.plot([x0, x1], [y0, y1], color=BAR_COLOR, linewidth=3.2, solid_capstyle='round', zorder=2, alpha=0.85)

# reference guide rings at 20/40/60% multicopy
REF_LABEL_ANGLE = math.radians(-63)  # empty whitespace in the lower-right, clear of every branch/label
for pct_ref in (20, 40, 60):
    r_ref = mc_r0 + RING_MAX_LEN * max_r * (pct_ref / MAX_PCT_FOR_SCALE)
    n_pts = 400
    xs = [r_ref * math.cos(2 * math.pi * i / n_pts) for i in range(n_pts + 1)]
    ys = [r_ref * math.sin(2 * math.pi * i / n_pts) for i in range(n_pts + 1)]
    ax.plot(xs, ys, color='#E4E8E5', linewidth=0.6, zorder=0, linestyle=':')
    lx, ly = r_ref * math.cos(REF_LABEL_ANGLE), r_ref * math.sin(REF_LABEL_ANGLE)
    ax.text(lx, ly, f'{pct_ref}%', fontsize=7, color='#9AA5A0', ha='center', va='bottom', zorder=0)

legend_handles = [
    mlines.Line2D([], [], color='#7A7A7A', linewidth=4.5, solid_capstyle='round',
                  label='Inner tick thickness: host-specific (single-phylum) cluster count at that order'),
    mlines.Line2D([], [], color=BAR_COLOR, linewidth=4.5, solid_capstyle='round',
                  label=f'Outer bar: % of that order’s phaC-positive genomes with >=2 copies (0-{MAX_PCT_FOR_SCALE:.0f}% scale)'),
]
ax.legend(handles=legend_handles, loc='upper left', fontsize=10, frameon=False)

# furthest actual content is the label anchor point (max_r + radius_step*offset);
# text glyphs extend further still depending on string length/rotation, so pad
# generously (same ~1.35x the base dendrogram itself uses for its own default limit).
label_r = max_r + RADIUS_STEP * LABEL_RADIUS_OFFSET
lim = label_r * 1.35
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

fig.suptitle('phaC Host-Specific Clusters by Order, with Multi-Copy Prevalence', fontsize=19, fontweight='bold', y=0.975)
fig.text(0.5, 0.935,
          f'{n_single:,} of {n_total:,} phaC clusters ({100*n_single/n_total:.1f}%) are found in exactly one host phylum — all {n_phyla} phyla shown, '
          f'top {TOP_K_CLASSES} classes/phylum, top {TOP_K_ORDERS} orders/class, {n_leaves} leaves total. Outer teal ring: genome-level multi-copy rate per leaf (section 9).',
          ha='center', fontsize=10.5, color='#5B6E70')
fig.text(0.5, 0.015,
          'Outer colored band = phylum. Inner tick thickness = host-specific cluster count at that order. Outer teal bars = % of that leaf’s own phaC-positive genomes with >=2 copies\n'
          '(same two QC corrections as section 9: no_hmm_triad_support exclusion + likely_artifact genomes dropped). Classes/orders beyond the top few per parent are folded into one\n'
          '"N smaller ..." leaf (grey) rather than omitted.',
          ha='center', fontsize=8.6, color='#5B6E70')

out_path = OUT / 'phac_taxonomy_radial_detailed_multicopy.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_taxonomy_radial_detailed_multicopy.pdf', facecolor='white')
