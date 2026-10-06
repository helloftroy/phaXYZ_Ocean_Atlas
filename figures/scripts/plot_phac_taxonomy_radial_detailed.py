"""Deeper version of plot_phac_taxonomy_radial.py: Domain -> Phylum ->
Class -> Order (two more layers than the original Phylum-only tree),
covering ALL phyla (no "N other phyla" bucket -- every phylum host-specific
(single-phylum) phaC clusters were found in gets its own branch, colored
ring wedge, and label).

Kingdom is NOT added above Domain: GTDB's real taxonomy has no rank above
domain (only Bacteria/Archaea/Eukarya exist as values), so a level there
would just be one trivial node splitting into <=3 children -- no genuine
structure. The actual fix for "just two branches with a bunch of twigs"
is adding depth belowClass, not width above Domain.

Density control: a flat "show every node above N clusters" threshold
doesn't control total leaf count well, because a handful of very
well-sampled phyla/classes (Pseudomonadota, Gammaproteobacteria...) can
still contribute hundreds of individually-qualifying orders. Instead this
caps each node's fan-out directly -- at most TOP_K_CLASSES per phylum and
TOP_K_ORDERS per shown class, the remainder folded into one "N smaller
..." leaf -- so total leaf count stays predictable regardless of how
taxonomically rich any single branch is.

QC: excludes phaC hits whose best_query is UNIPROT:C7BNH2, a mislabeled
reference (actually isochorismate synthase, not a PHA synthase) -- see
_phac_qc.py. A cluster whose entire membership was C7BNH2-only hits drops
out of the tree entirely (correctly -- it was never really a phaC cluster).
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _radial_dendrogram as rd
import _phac_qc
import _phylum_colors

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Wedge
import colorsys

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
OUT = Path(__file__).resolve().parent.parent
TOP_K_CLASSES = 3
TOP_K_ORDERS = 2
OTHER_COLOR = '#7A7A7A'

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

# per-phylum -> per-class -> per-order counts (all real, un-thresholded)
phylum_totals = defaultdict(int)
phylum_domain = {}
class_counts_by_phylum = defaultdict(lambda: defaultdict(int))
order_counts_by_class = defaultdict(lambda: defaultdict(int))  # (phylum,class) -> order -> count

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
            paths.append(((dom, phylum, cls, order), order_n))
        if other_orders:
            n_other = len(other_orders)
            other_n = sum(n for _, n in other_orders)
            paths.append(((dom, phylum, cls, f'{n_other} smaller orders'), other_n))

    if other_classes:
        n_other_cls = len(other_classes)
        other_cls_n = sum(n for _, n in other_classes)
        paths.append(((dom, phylum, f'{n_other_cls} smaller classes', '(mixed)'), other_cls_n))

root = rd.build_tree(paths)
n_leaves = len(paths)
print(f'{n_leaves} leaves shown (order-level, or "smaller classes/orders" buckets)')


# Phylum colours come from _phylum_colors so the tree, the prevalence bubble
# plot and the synthase-class bars agree; they are assembled into one figure.
# The generator that used to sit here moved there unchanged, so the tree looks
# the same as before.
phyla_ranked = sorted(phylum_totals.items(), key=lambda x: -x[1])
PHYLUM_COLORS = {p: _phylum_colors.color_for(p) for p, _ in phyla_ranked}

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 11})
fig, ax = plt.subplots(figsize=(12, 12), dpi=300)
fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

def branch_color(path, count):
    # raises IndexError for depth-1 (domain) paths -- caught by the
    # renderer, which falls back to the flat default edge_color for those,
    # since a domain spans many differently-colored phyla below it.
    if 'smaller' in path[-1] or 'smaller' in path[-2]:
        return OTHER_COLOR
    return PHYLUM_COLORS.get(path[1], OTHER_COLOR)


# Compress the radial depth with non-uniform radii: keep the terminal leaves
# and phylum ring readable, but pull the internal taxonomic levels inward so
# the center does not balloon into a sparse, long-limbed circle.
# Arms shortened again for the manuscript panel: the terminal ring moves in from
# 0.40 to 0.33 and the internal levels compress with it, so the leaf spokes stop
# dominating and the phylum labels can grow without the figure getting bigger.
RADII_BY_DEPTH = [0.0, 0.045, 0.105, 0.185, 0.300]
RADIUS_STEP = 0.56
LABEL_RADIUS_OFFSET = 0.04
LEAF_LABEL_MIN_COUNT = 10**9


def leaf_label(path, count):
    label = path[-1]
    if count >= LEAF_LABEL_MIN_COUNT and 'smaller' not in label and label != '(mixed)':
        return label if len(label) <= 22 else label[:20] + '…'
    return ''

angle_of, max_r = rd.draw_radial_dendrogram(
    ax, root, radius_step=RADIUS_STEP, edge_color='#D8DED7', edge_width=0.7,
    leaf_label_fn=leaf_label, label_levels=(4,), label_fontsize=5.5, label_radius_offset=LABEL_RADIUS_OFFSET,
    leaf_size_fn=lambda path, count: 6 + 70 * (count / max(n for _, n in paths)),
    leaf_color_fn=branch_color, edge_color_fn=branch_color,
    show_internal_dots=True, internal_dot_size=4.2, internal_dot_color='#DDE2DB',
    radii_by_depth=RADII_BY_DEPTH,
)

# ---- outer colored clade ring, one wedge per phylum (depth 2) ----
# Keep the phylum band close to the terminal leaf radius. The previous pass
# shortened the tree but left this band much farther out, creating a large
# empty annulus between leaf dots and taxonomy blocks.
ring_r0 = 0.322
ring_r1 = 0.372
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
    if span_deg > 1.3:  # only label wedges wide enough to hold readable text
        mid = (lo + hi) / 2
        lx, ly = rd._polar(mid, ring_r1 + 0.028)
        rot = math.degrees(mid)
        ha = 'left'
        if math.cos(mid) < 0:
            rot += 180
            ha = 'right'
        short_label = label if len(label) < 20 else label[:18] + '…'
        fontsize = 13.5 if span_deg > 6 else 9.5
        ax.text(lx, ly, short_label, rotation=rot, rotation_mode='anchor', ha=ha, va='center',
                fontsize=fontsize, fontweight='bold', color='#20302C', zorder=5,
                path_effects=[pe.withStroke(linewidth=3.0, foreground='white')])

lim = 0.78
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

# No title, subtitle or footnote: this panel goes into a manuscript figure.
# bbox_inches='tight' because the drawn content is a circle inscribed in a square
# axis: shortening the arms shrank the circle but left the canvas the same size,
# so the panel arrived with a wide white border to crop by hand.
out_path = OUT / 'phac_taxonomy_radial_detailed.png'
fig.savefig(out_path, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.04)
print('saved', out_path)
fig.savefig(OUT / 'phac_taxonomy_radial_detailed.pdf', facecolor='white', bbox_inches='tight', pad_inches=0.04)
