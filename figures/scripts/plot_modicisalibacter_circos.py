"""Chord diagram (circos-style, not a genome-circular plot): one sector
per genome, that genome's own phaC copies placed as nodes around its
sector, and chords connecting phaC copies from DIFFERENT genomes that
belong to the same 70%-identity paralog cluster -- i.e. "this genome's
copy of paralog X connects to that genome's copy of the same paralog X".
Chord opacity/width scales with the exact pairwise %identity between the
two specific proteins (computed directly, not read off the cluster
threshold); node color encodes which paralog cluster a protein belongs to
(the same 70%-cluster framework used throughout section 9). Prototype
case, picked deliberately small and already well-understood: Modicisalibacter
zincidurans (section 9.5's 6 genomes) plus 4 "outside" genomes -- one
representative each from Cobetia, Vreelandella, Halomonas, and
Marinobacter, the four genera section 9.5's addendum found sharing
Modicisalibacter's paralog clusters. Same method is meant to extend to
HK1/CARD22-1 next (piece 5 of the current request) once this is validated
on a small, checkable case.

Outside-genome selection: for each of the two Modicisalibacter clusters
with real cross-genus membership (...602748, Halomonadaceae-specific;
...783925, the broader cross-family lineage), the highest-completeness/
lowest-contamination genome of each contributing genus was picked (ties
broken toward completeness) -- KOPF15-1_SAMEA7392432_MAG_00000114
(Cobetia), OBRI24-1_SAMEA115094377_MAG_00000066 (Vreelandella),
RSGB23-1_GCF-003547075-V1_GENO_10000001 (Halomonas, an NCBI isolate),
SCHR19-1_SAMN09405896_MAG_00000020 (Marinobacter). Each contributes its
FULL phaC complement to the diagram, not just the one target that placed
it in a shared cluster -- so their own additional, Modicisalibacter-
unconnected paralogs show up too (5 more distinct clusters appear this
way, unconnected singletons in the diagram, which is itself informative:
Modicisalibacter's own paralog repertoire is not simply a subset of what
these other genera carry).

Usage:
    python figures/scripts/plot_modicisalibacter_circos.py

Outputs:
    figures/modicisalibacter_circos.png / .pdf
    figures/modicisalibacter_circos_nodes.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

MODZ_GENOMES = [
    'DUAR15-1_SAMN06266142_MAG_00000002', 'RSGB23-1_GCF-000731955-V1_GENO_10000001',
    'SANC23-1_SAMEA110646468_MAG_00000065', 'TARA_SAMEA2623562_MAG_00000025',
    'TARA_SAMEA2623564_MAG_00000048', 'ZHEN20-1_SAMN07748058_MAG_00000118',
]
OUTSIDE_GENOMES = {
    'Cobetia': 'KOPF15-1_SAMEA7392432_MAG_00000114',
    'Vreelandella': 'OBRI24-1_SAMEA115094377_MAG_00000066',
    'Halomonas': 'RSGB23-1_GCF-003547075-V1_GENO_10000001',
    'Marinobacter': 'SCHR19-1_SAMN09405896_MAG_00000020',
}
GENOME_ORDER = MODZ_GENOMES + list(OUTSIDE_GENOMES.values())
GENOME_LABEL = {g: f'M.z. {g.split("_")[0]}-{g.split("_")[1][-4:]}' for g in MODZ_GENOMES}
GENOME_LABEL.update({g: f'{genus} ({g.split("_")[0]})' for genus, g in OUTSIDE_GENOMES.items()})
GENOME_IS_OUTSIDE = {g: (g in OUTSIDE_GENOMES.values()) for g in GENOME_ORDER}

bad_targets = _phac_qc.load_bad_targets()

# ---------------------------------------------------------------------
# 1. per-genome target lists (corrected NR100 join)
# ---------------------------------------------------------------------
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in GENOME_ORDER and row['target_id'] not in bad_targets:
            genome_targets[row['genome']].add(row['target_id'])

all_tids = set()
for s in genome_targets.values():
    all_tids |= s
print(f'{len(GENOME_ORDER)} genomes, {sum(len(v) for v in genome_targets.values())} protein-genome instances, '
      f'{len(all_tids)} distinct target_ids')

# ---------------------------------------------------------------------
# 2. cluster assignment, sequences, protein length
# ---------------------------------------------------------------------
member_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if row[1] in all_tids:
            member_to_cluster[row[1]] = row[0]

seqs = {}
cur_id, cur_seq = None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in all_tids:
                seqs[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in all_tids:
        seqs[cur_id] = ''.join(cur_seq)

# ---------------------------------------------------------------------
# 3. build node list: one node per (genome, target_id) instance
# ---------------------------------------------------------------------
nodes = []  # dicts: genome, target_id, cluster, angle_index_within_sector
for g in GENOME_ORDER:
    for t in sorted(genome_targets[g]):
        nodes.append({'genome': g, 'target_id': t, 'cluster': member_to_cluster.get(t, t)})

# color palette per distinct cluster (paralog group), consistent across
# both Modicisalibacter's own 4 known clusters and any new ones the
# outside genomes bring in
cluster_ids = sorted(set(n['cluster'] for n in nodes), key=lambda c: -sum(1 for n in nodes if n['cluster'] == c))
PALETTE = ['#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B', '#7A5FA0', '#C9A227', '#4A7FB5',
           '#B33951', '#6FA88A', '#8B4513', '#2E86AB', '#D98E04']
cluster_color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(cluster_ids)}
KNOWN_MODZ_CLUSTER_NAME = {
    'OMDBv2.0_AA_G_NR100_000148602748': 'paralog 602748 (Class I)',
    'OMDBv2.0_AA_G_NR100_000139783925': 'paralog 783925 (phaC/phaZ-ambig.)',
    'OMDBv2.0_AA_G_NR100_000061282509': 'paralog 282509 (98.8% struct. id.)',
    'OMDBv2.0_AA_G_NR100_000139753017': 'paralog 753017 (738aa, atypical)',
}


def cluster_label(c):
    return KNOWN_MODZ_CLUSTER_NAME.get(c, f'other ({c[-6:]})')


# ---------------------------------------------------------------------
# 4. pairwise %identity for every distinct sequence pair (chord weight)
# ---------------------------------------------------------------------
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = 'global'


def pident(a, b):
    aln = aligner.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    m = sum(1 for x, y in zip(s1, s2) if x == y and x != '-')
    al = sum(1 for x, y in zip(s1, s2) if x != '-' and y != '-')
    return 100 * m / al if al else 0.0


distinct_tids = sorted(seqs)
pident_cache = {}
for i, t1 in enumerate(distinct_tids):
    for t2 in distinct_tids[i + 1:]:
        pident_cache[frozenset((t1, t2))] = pident(seqs[t1], seqs[t2])
print(f'{len(pident_cache)} pairwise identities computed')

# ---------------------------------------------------------------------
# 5. write node table
# ---------------------------------------------------------------------
nodes_out = OUT / 'modicisalibacter_circos_nodes.tsv'
with open(nodes_out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['genome', 'genome_label', 'target_id', 'cluster', 'cluster_label', 'protein_length_aa', 'is_outside_genome'])
    for n in nodes:
        w.writerow([n['genome'], GENOME_LABEL[n['genome']], n['target_id'], n['cluster'], cluster_label(n['cluster']),
                    len(seqs.get(n['target_id'], '')), GENOME_IS_OUTSIDE[n['genome']]])
print(f'wrote {nodes_out}')

# ---------------------------------------------------------------------
# plotting -- circle geometry
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
fig, ax = plt.subplots(figsize=(13, 13), dpi=300, subplot_kw={'aspect': 'equal'})

GAP_DEG = 3.0   # gap between genome sectors
R_OUTER = 1.0   # sector arc radius
R_NODE = 0.93   # node marker radius
R_LABEL = 1.08  # label radius
R_CHORD_CTRL = 0.0  # bezier control point at center -> classic chord-diagram curve

n_genomes = len(GENOME_ORDER)
total_gap = GAP_DEG * n_genomes
counts = {g: max(len(genome_targets[g]), 1) for g in GENOME_ORDER}
total_count = sum(counts.values())
avail_deg = 360.0 - total_gap

sector_span = {}
start = 90.0  # start at top, go clockwise
sector_start = {}
for g in GENOME_ORDER:
    span = avail_deg * counts[g] / total_count
    sector_start[g] = start
    sector_span[g] = span
    start -= (span + GAP_DEG)

node_angle = {}  # (genome, target_id) -> angle in degrees
for g in GENOME_ORDER:
    tids_g = sorted(genome_targets[g])
    span = sector_span[g]
    n = len(tids_g)
    for i, t in enumerate(tids_g):
        # evenly spaced within the sector, centered
        frac = (i + 0.5) / n
        node_angle[(g, t)] = sector_start[g] - frac * span


def polar(r, deg):
    rad = np.radians(deg)
    return r * np.cos(rad), r * np.sin(rad)


# sector arcs (background bands), colored by Modicisalibacter vs outside genus
GENUS_BAND_COLOR = {**{g: '#DCEDEA' for g in MODZ_GENOMES},
                     OUTSIDE_GENOMES['Cobetia']: '#F3E3D3', OUTSIDE_GENOMES['Vreelandella']: '#EFD9E8',
                     OUTSIDE_GENOMES['Halomonas']: '#E3E9F7', OUTSIDE_GENOMES['Marinobacter']: '#EAF0DC'}
for g in GENOME_ORDER:
    wedge = mpatches.Wedge((0, 0), R_OUTER + 0.05, sector_start[g] - sector_span[g], sector_start[g],
                             width=0.09, facecolor=GENUS_BAND_COLOR.get(g, '#E5E5E5'), edgecolor='#8B8F8C', linewidth=0.6, zorder=2)
    ax.add_patch(wedge)
    mid = sector_start[g] - sector_span[g] / 2
    x, y = polar(R_LABEL, mid)
    mid_norm = (mid + 180) % 360 - 180  # normalize to (-180, 180] for the half-check below
    rot = mid_norm if -90 < mid_norm <= 90 else mid_norm + 180
    ha = 'left' if -90 < mid_norm <= 90 else 'right'
    fw = 'bold' if GENOME_IS_OUTSIDE[g] else 'normal'
    ax.text(x, y, GENOME_LABEL[g], rotation=rot, ha=ha, va='center', fontsize=8.3, fontweight=fw,
             rotation_mode='anchor', color='#20302C')

# chords: connect same-cluster nodes across DIFFERENT genomes
drawn = set()
n_chords = 0
for i, n1 in enumerate(nodes):
    for n2 in nodes[i + 1:]:
        if n1['genome'] == n2['genome']:
            continue
        if n1['cluster'] != n2['cluster']:
            continue
        key = frozenset((n1['target_id'], n2['target_id'], n1['genome'], n2['genome']))
        if key in drawn:
            continue
        drawn.add(key)
        pid = 100.0 if n1['target_id'] == n2['target_id'] else pident_cache.get(frozenset((n1['target_id'], n2['target_id'])), 0)
        a1 = node_angle[(n1['genome'], n1['target_id'])]
        a2 = node_angle[(n2['genome'], n2['target_id'])]
        x1, y1 = polar(R_NODE, a1)
        x2, y2 = polar(R_NODE, a2)
        path = MplPath([(x1, y1), (0, 0), (x2, y2)], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
        alpha = 0.15 + 0.75 * (pid / 100)
        lw = 0.5 + 3.0 * (pid / 100)
        color = cluster_color[n1['cluster']]
        patch = mpatches.PathPatch(path, facecolor='none', edgecolor=color, linewidth=lw, alpha=alpha, zorder=1)
        ax.add_patch(patch)
        n_chords += 1
print(f'{n_chords} chords drawn')

# node markers
for n in nodes:
    a = node_angle[(n['genome'], n['target_id'])]
    x, y = polar(R_NODE, a)
    ax.scatter([x], [y], s=70, color=cluster_color[n['cluster']], edgecolor='#20302C', linewidth=0.8, zorder=3)

ax.set_xlim(-1.35, 1.35)
ax.set_ylim(-1.35, 1.35)
ax.axis('off')

legend_handles = [mpatches.Patch(color=cluster_color[c], label=cluster_label(c)) for c in cluster_ids
                   if sum(1 for n in nodes if n['cluster'] == c) > 1 or c in KNOWN_MODZ_CLUSTER_NAME]
legend_handles.append(mpatches.Patch(color='#B0B6B2', label='singleton (no shared paralog cluster in this set)'))
legend1 = ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.28, -0.04), fontsize=8, frameon=False, ncol=1,
           title='Paralog group (node color)', title_fontsize=8.5)
ax.add_artist(legend1)  # preserve this legend before the second ax.legend() call below replaces it

band_handles = [mpatches.Patch(color='#DCEDEA', label='Modicisalibacter zincidurans (6 genomes)')]
for genus, color in [('Cobetia', '#F3E3D3'), ('Vreelandella', '#EFD9E8'), ('Halomonas', '#E3E9F7'), ('Marinobacter', '#EAF0DC')]:
    band_handles.append(mpatches.Patch(color=color, label=f'{genus} (outside genus)'))
ax.legend(handles=band_handles, loc='upper center', bbox_to_anchor=(0.72, -0.04), fontsize=8, frameon=False,
           title='Genome sector', title_fontsize=8.5)

fig.suptitle('phaC paralogs across Modicisalibacter zincidurans and its closest outside relatives',
             fontsize=15, fontweight='bold', y=0.98)
fig.text(0.5, 0.945, 'Each sector = one genome; each dot = one phaC copy. Chords connect copies from different genomes '
          'sharing the same 70%-identity paralog cluster;\nchord opacity/width = exact pairwise %identity between those two proteins.',
          ha='center', fontsize=9, color='#5B6E70')

out_path = OUT / 'modicisalibacter_circos.png'
fig.savefig(out_path, dpi=300, facecolor='white', bbox_inches='tight', pad_inches=0.35)
print('\nsaved', out_path)
fig.savefig(OUT / 'modicisalibacter_circos.pdf', facecolor='white', bbox_inches='tight', pad_inches=0.35)
print('saved pdf too')
