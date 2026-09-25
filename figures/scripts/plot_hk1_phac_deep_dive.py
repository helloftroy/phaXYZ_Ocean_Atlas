"""Unknown-HK1 deep dive -- section 9.5 picked Modicisalibacter zincidurans
over this lineage as the biotech-potential case study, but never actually
looked closely at HK1's own single highest-copy genome the way section
9.4 (CARD22-1) and 9.5 (Modicisalibacter) did. Picked up directly on
request, scoped to just two things (per that request -- the NJ-tree/
structural-class-HMM parts of the older deep dives are now superseded by
real ESMFold structures, section 9.11-9.12): a %identity heatmap and a
contig map.

Genome: CHAS20-1_SAMN14350117_MAG_00000061, the single highest-copy HK1
genome (7 phaC copies, matching section 9.5's own "had a 7-copy genome"
mention of this lineage) -- Pacific marine sediment (-17.47, -149.78),
77.5% complete / 8.5% contamination (below the 10% likely_artifact gate,
section 9.3), n_families_present=10 (phaA=9, phaB=16, phaC=7, phaJ=7,
phaZ=9, phaG/phaY/phaP/phaF=1 each, phaR_regulator=2, phaE=0 -- a rich
PHA gene repertoire, but no phaE, ruling out Class III's phaE-dependent
heterodimer mechanism for this genome specifically).

All 7 copies have a folded PDB (unlike CARD22-1/Modicisalibacter's
partial coverage) and all 7 are structurally triad-complete by
find_structural_triad.py's geometric check -- including target ...711722
(737aa, the longest/most divergent copy, matches the same weakly-
annotated Q5P962 reference flagged in section 9.7/9.12), which the older
alignment-column method calls triad_complete=False -- the same
"hmmalign misses it in long divergent sequences" pattern already found
in CARD22-1's 708aa copy. Structurally, all 7 look like real, intact
active sites.

Contig map caveat, stated directly rather than silently reusing section
9.4/9.5's method: those two deep dives' "other pha genes on this
genome's scaffolds" panels required a one-off fetch of this genome's own
prodigal gene calls from the OMDB mirror (individual per-genome
.genes.faa.gz files). Checked directly for this genome (2026-09-25): the
data provider has since reorganized that mirror (as of the directory's
own 2026-08-25 timestamp) -- individual per-genome gene-call files are
no longer served; only one ~4.6GB .tar archive per whole study (here,
CHAS20-1.tar) is available, impractical to download just for one
genome's gene positions. This map therefore shows each phaC copy's own
scaffold + gene index (available locally, via phaC_cluster_sequences.faa's
"rep=" field) without the neighboring-gene overlay the CARD22-1/
Modicisalibacter maps have. That absence is itself informative here: all
7 copies sit on 7 DIFFERENT scaffolds (no two share one), so even with
neighbor annotation this genome was never going to show the kind of
tight phaE-phaC-phaJ or phaB-phaA-phaC operon those two other genomes did
-- a fragmented-assembly-consistent, one-gene-per-scaffold picture is
already the visible result of the data actually available.

Usage:
    python figures/scripts/plot_hk1_phac_deep_dive.py

Outputs:
    figures/hk1_phac_identity_heatmap.png / .pdf
    figures/hk1_phac_contig_map.png / .pdf
    figures/hk1_phac_paralogs.tsv
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / 'catalytic_domain'))
import _phac_qc
from find_structural_triad import structural_triad_complete_batch

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

GENOME = 'CHAS20-1_SAMN14350117_MAG_00000061'

# ---------------------------------------------------------------------
# 1. this genome's own phaC targets, current QC applied
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
targets = []
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] == GENOME and row['target_id'] not in bad_targets:
            targets.append(row['target_id'])
print(f'{GENOME}: {len(targets)} QC-passing phaC targets')

# own_scaffold/own_gene_index -- from phaC_cluster_sequences.faa's "rep=" field,
# restricted to entries where the rep IS this genome (see module docstring
# for why this covers phaC's own position but not neighboring gene families)
own_pos = {}
seqs = {}
wanted = set(targets)
cur_id, cur_seq = None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in wanted:
                seqs[cur_id] = ''.join(cur_seq)
            parts = line[1:].split()
            cur_id = parts[0]
            cur_seq = []
            if len(parts) > 1 and parts[1].startswith('rep=') and cur_id in wanted:
                rep = parts[1][len('rep='):].split(';')[0]
                if rep.startswith(GENOME + '-'):
                    scaf_gene = rep[len(GENOME) + 1:]
                    scaf, gene_idx = scaf_gene.rsplit('_', 1)
                    own_pos[cur_id] = (scaf, int(gene_idx))
        else:
            cur_seq.append(line.strip())
    if cur_id in wanted:
        seqs[cur_id] = ''.join(cur_seq)

triad_complete_alignment = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in wanted:
            triad_complete_alignment[row['target_id']] = row['triad_complete'] == 'True'

triad_complete_structural = structural_triad_complete_batch(targets)

meta = {}
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in wanted:
            meta[row['target_id']] = row

targets.sort(key=lambda t: own_pos.get(t, ('', 0)))

rows_out = []
for t in targets:
    scaf, idx = own_pos.get(t, ('unknown', -1))
    rows_out.append({
        'target_id': t, 'protein_length_aa': len(seqs.get(t, '')),
        'own_scaffold': scaf, 'own_gene_index': idx,
        'best_query': meta.get(t, {}).get('best_query', ''),
        'pident_to_best_ref': meta.get(t, {}).get('best_pident', ''),
        'triad_complete_structural': triad_complete_structural.get(t),
        'triad_complete_alignment_column': triad_complete_alignment.get(t),
    })
paralogs_out = OUT / 'hk1_phac_paralogs.tsv'
with open(paralogs_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)
print('wrote', paralogs_out)

# ---------------------------------------------------------------------
# 2. %identity heatmap (BLOSUM62 global pairwise, same method used
#    throughout section 9)
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


n = len(targets)
mat = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        mat[i, j] = 100.0 if i == j else pident(seqs[targets[i]], seqs[targets[j]])
print(f'{n*(n-1)//2} pairwise identities computed')

short_labels = [f'{t[-6:]}\n({len(seqs[t])}aa)' for t in targets]

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
fig, ax = plt.subplots(figsize=(8.5, 7.5), dpi=300)
im = ax.imshow(mat, cmap='YlGnBu', vmin=0, vmax=100)
ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(short_labels, fontsize=8.5, rotation=45, ha='right')
ax.set_yticklabels(short_labels, fontsize=8.5)
for i in range(n):
    for j in range(n):
        val = mat[i, j]
        color = 'white' if val > 60 else '#20302C'
        ax.text(j, i, f'{val:.0f}', ha='center', va='center', fontsize=8.5, color=color)
ax.set_title(f'phaC pairwise %identity: {GENOME}\n(7 copies, "Unknown HK1")', fontsize=12.5, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('%identity (BLOSUM62 global alignment)')
fig.text(0.5, 0.01, 'All 7 copies are structurally triad-complete (find_structural_triad.py, section 9.12) despite low pairwise identity to each other --\n'
                     'independent acquisitions/divergence, not a single recently-duplicated family, the same pattern section 9.4 found in CARD22-1.',
          ha='center', fontsize=8, color='#5B6E70')
fig.tight_layout(rect=[0, 0.045, 1, 1])
fig.savefig(OUT / 'hk1_phac_identity_heatmap.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'hk1_phac_identity_heatmap.pdf', facecolor='white')
print('saved', OUT / 'hk1_phac_identity_heatmap.png')

# ---------------------------------------------------------------------
# 3. contig map -- one short track per scaffold (all 7 copies are on
#    7 DIFFERENT scaffolds -- see module docstring)
# ---------------------------------------------------------------------
INK = '#20302C'
MUTED = '#667879'
PHAC_COLOR = '#2F6F63'

fig2, axes = plt.subplots(n, 1, figsize=(9.5, 1.0 * n + 0.9), dpi=300)
for ax_i, t in zip(axes, targets):
    scaf, idx = own_pos[t]
    ax_i.set_xlim(-1, 1)
    ax_i.set_ylim(-0.65, 0.7)
    ax_i.axis('off')
    ax_i.plot([-0.7, 0.7], [0, 0], color='#B9C7C1', linewidth=1.2, zorder=1)
    arrow = FancyArrowPatch((-0.16, 0), (0.16, 0),
                              arrowstyle='Simple,tail_width=0.5,head_width=0.95,head_length=0.4',
                              mutation_scale=15, linewidth=1.1, facecolor=PHAC_COLOR, edgecolor=INK, zorder=3)
    ax_i.add_patch(arrow)
    triad_s = triad_complete_structural.get(t)
    triad_tag = 'triad-complete (structural)' if triad_s else ('triad-complete? (alignment-column only, no PDB)' if triad_s is None else 'not triad-complete')
    ax_i.text(0, 0.32, f'phaC {t[-6:]}  --  {len(seqs[t])}aa', ha='center', va='bottom', fontsize=9, fontweight='bold', color=INK)
    ax_i.text(0, -0.30, triad_tag, ha='center', va='top', fontsize=7.3, color=MUTED)
    ax_i.text(-0.95, 0, f'{scaf}\ngene {idx}', ha='left', va='center', fontsize=8, color=INK, fontweight='bold')

fig2.suptitle(f'phaC contig positions: {GENOME}', fontsize=13, fontweight='bold', y=0.99)
fig2.text(0.5, 0.945, 'All 7 copies sit on 7 different scaffolds -- no operon-style co-location (unlike CARD22-1/Modicisalibacter)',
          ha='center', fontsize=9.5, color='#3A4A46')
fig2.text(0.5, 0.008, 'Neighboring phaA/B/E/J/Z genes not shown: the OMDB mirror no longer serves individual per-genome gene-call files\n'
                       '(reorganized into ~GB-scale per-study .tar archives since the CARD22-1/Modicisalibacter maps were built) -- see script docstring.',
          ha='center', fontsize=7.6, color='#5B6E70')
fig2.subplots_adjust(left=0.17, right=0.97, top=0.85, bottom=0.06, hspace=0.35)
fig2.savefig(OUT / 'hk1_phac_contig_map.png', dpi=300, facecolor='white')
fig2.savefig(OUT / 'hk1_phac_contig_map.pdf', facecolor='white')
print('saved', OUT / 'hk1_phac_contig_map.png')
