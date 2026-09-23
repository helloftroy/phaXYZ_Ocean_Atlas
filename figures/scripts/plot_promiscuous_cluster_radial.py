"""Radial dendrogram of host taxonomic spread for ONE phaC cluster --
built for the two most taxonomically promiscuous clusters, re-derived
twice now: once after the original reference-query audit (67 accessions),
and again after the deeper 2026-09-22 phaC fix (219 accessions total --
see figures/PHA_CLEAN_RESULTS.md section 2). Both earlier "most
promiscuous" picks have now disappeared entirely at least once each:
...131018176 and ...095227834 (the pre-audit picks) collapsed to
near-nothing after the first audit; ...246448549 (REDSEA-S09-B13, the
post-first-audit #1 pick) now has ZERO genomes after the second fix --
its "most promiscuous cluster" status turns out to have been driven
entirely by targets recruited from among the 152 newly-excluded
references, not a real biological pattern at any point. The current
top-2 by distinct-phyla count (computed fresh from the corrected
cluster/metadata files) are ...117564048 (713 genomes, 28 genera, 8
phyla -- archaea-dominated, Nitrosopumilus-heavy) and ...022439766 (380
genomes, 27 genera, 8 phyla -- also archaea-dominated, Nitrosopelagicus-
heavy; genome/genera counts shifted from the prior pass but this cluster
has now survived BOTH rounds of correction). Genus is rendered as
unlabeled leaf dots (too many to label legibly at this scale) sized by
genome count within THIS cluster; phylum is the labeled outer ring --
same design as plot_phac_taxonomy_radial_detailed.py (colored ring,
branches colored to match the leaves they lead to, labels pushed past
the ring rather than run through it).

QC: excludes every reference confirmed wrong-gene in the full audit
(figures/PHA_CLEAN_RESULTS.md section 2) -- see _phac_qc.py.
"""
import colorsys
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
from matplotlib.patches import Wedge

FA = Path('/Users/hellpark/multimodal_seusmbol/PHA_Ocean_Atlas/PHA_bioprospecting/omdb_search/results')
OUT = Path(__file__).resolve().parent.parent
OTHER_COLOR = '#7A7A7A'
RADIUS_STEP = 1.15

CLUSTERS = [
    ('OMDBv2.0_AA_G_NR100_000117564048', 'Nitrosopumilus'),
    ('OMDBv2.0_AA_G_NR100_000022439766', 'Nitrosopelagicus'),
]

assignments = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        assignments[row[1]] = row[0]


def gen_phylum_colors(phyla_sorted_by_size):
    """Golden-angle hue stepping -- see plot_phac_taxonomy_radial_detailed.py
    for why (evenly-spaced-by-rank hues look like a gradient, not a real
    categorical palette, since rank-adjacent phyla are usually also
    spatially adjacent in the tree)."""
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


for full_cid, dominant_genus in CLUSTERS:
    short_id = full_cid[-12:]
    genus_counts = defaultdict(int)
    genus_phylum = {}
    genus_domain = {}
    seen_genomes = set()

    with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if _phac_qc.is_bad(row['best_query']):
                continue
            if assignments.get(row['target_id']) != full_cid:
                continue
            genome = row['genome']
            if genome in seen_genomes:
                continue
            seen_genomes.add(genome)
            genus = row.get('gtdb_genus', '') or 'Unknown'
            genus_counts[genus] += 1
            genus_phylum[genus] = row.get('gtdb_phylum', '') or 'Unknown'
            genus_domain[genus] = row.get('gtdb_domain', '') or 'Bacteria'

    n_genomes = len(seen_genomes)
    n_genera = len(genus_counts)
    phylum_totals = defaultdict(int)
    for g, n in genus_counts.items():
        phylum_totals[genus_phylum[g]] += n
    n_phyla = len(phylum_totals)
    print(f'{short_id}: {n_genomes} genomes, {n_genera} genera, {n_phyla} phyla')

    phyla_ranked = sorted(phylum_totals.items(), key=lambda x: -x[1])
    PHYLUM_COLORS = gen_phylum_colors([p for p, _ in phyla_ranked])

    def branch_color(path, count):
        # raises IndexError for depth-1 (domain) paths -- caught by the
        # renderer, falling back to the flat default edge_color, since a
        # domain spans multiple differently-colored phyla below it.
        return PHYLUM_COLORS.get(path[1], OTHER_COLOR)

    paths = [((genus_domain[g], genus_phylum[g], g), n) for g, n in genus_counts.items()]
    root = rd.build_tree(paths)

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, ax = plt.subplots(figsize=(14, 14), dpi=300)

    # leaf (genus) labels are never drawn (too many); label_radius_offset
    # is irrelevant here since label_levels is empty, but kept consistent
    # with the detailed-taxonomy script's approach for the ring math below
    angle_of, max_r = rd.draw_radial_dendrogram(
        ax, root, radius_step=RADIUS_STEP, edge_color='#D8DED7', edge_width=0.8,
        label_levels=(), leaf_size_fn=lambda path, count: 14 + 90 * (count / max(genus_counts.values())),
        leaf_color_fn=branch_color, edge_color_fn=branch_color,
    )

    # ---- outer colored clade ring, one wedge per phylum (depth 2) ----
    ring_r0 = max_r * 1.05
    ring_r1 = max_r * 1.13
    phylum_nodes = []
    for dom_label, dom_child in rd._sorted_children(root[1]):
        for label, child in rd._sorted_children(dom_child[1]):
            phylum_nodes.append((label, child))

    for label, child in phylum_nodes:
        lo, hi = rd.leaf_angle_range(child, angle_of)
        pad = 0.004
        color = PHYLUM_COLORS.get(label, OTHER_COLOR)
        wedge = Wedge((0, 0), ring_r1, math.degrees(lo - pad), math.degrees(hi + pad),
                      width=ring_r1 - ring_r0, facecolor=color, alpha=0.85, edgecolor='white', linewidth=0.8, zorder=0)
        ax.add_patch(wedge)
        mid = (lo + hi) / 2
        lx, ly = rd._polar(mid, ring_r1 + max_r * 0.06)
        rot = math.degrees(mid)
        ha = 'left'
        if math.cos(mid) < 0:
            rot += 180
            ha = 'right'
        ax.text(lx, ly, label, rotation=rot, rotation_mode='anchor', ha=ha, va='center',
                fontsize=11, fontweight='bold', color='#20302C', zorder=5,
                path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])

    lim = ring_r1 * 1.3
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    fig.suptitle(f'phaC Cluster …{short_id} ({dominant_genus}-dominant): Host Spread', fontsize=17, fontweight='bold', y=0.975)
    fig.text(0.5, 0.93, f'{n_genomes:,} genomes across {n_genera} genera and {n_phyla} phyla — one dot per genus, sized by genome count',
              ha='center', fontsize=10.5, color='#5B6E70')
    fig.text(0.5, 0.015, 'Outer colored band = phylum, matching branch/dot color. Genus-level dots are intentionally unlabeled (too many to read individually) -- '
              'this is one of the two most taxonomically promiscuous clusters in the dataset.',
              ha='center', fontsize=8.6, color='#5B6E70')

    out_path = OUT / f'promiscuous_cluster_{short_id}_radial.png'
    fig.savefig(out_path, dpi=300, facecolor='white')
    print('saved', out_path)
    fig.savefig(OUT / f'promiscuous_cluster_{short_id}_radial.pdf', facecolor='white')
    plt.close(fig)
