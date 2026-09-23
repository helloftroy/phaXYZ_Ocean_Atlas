"""Single-genome deep dive, picked from section 9.3's scripted audit
survivors: CARD22-1_SAMN24292811_MAG_00000010, an uncultured Desulfoluna
MAG from Red Sea coral tissue (Goniastrea edwardsi, 22.18N 39.02E, study
CARD22-1) -- 7 phaC copies, all 7 in distinct 70%-identity clusters,
99.35% CheckM-complete / 1.77% contamination, `likely_legit` call. Chosen
over the higher-raw-copy-count RSGB23-1/Legionella cases specifically
because it is a real environmental sample at a specific location with a
clean bin, not a location-blind isolate batch (section 9.2) or the kind
of fragmented assembly section 9.2 flagged as an artifact risk.

Four questions, each answered with a different local data source:

1. Do the copies form one clade (recent duplication) or fall into
   different classes (independent acquisitions)? -- NJ tree + pairwise
   %identity, computed directly from the 7 sequences with Biopython's
   pairwise aligner (BLOSUM62, global) and Bio.Phylo's distance-tree
   constructor -- no external aligner needed, same method as section
   9.2's Legionella tree.

2. Contig map, phaC and the OTHER pha genes together -- the corrected
   NR100 target_id<->genome join gives phaC scaffold positions directly
   (via phaC_cluster_sequences.faa's "rep=" field, restricted to entries
   where the rep IS this genome). For the other pha families (phaA/B/E/
   J/Z), no equivalent position-annotated per-family sequence catalog
   exists in this repo (only phaC has one) and this genome's raw prodigal
   gene calls are not cached locally, so OTHER_PHA_GENES_ON_SCAFFOLDS
   below is a one-off, manually-run, fully-documented lookup: (a) fetched
   this genome's own gene calls directly from the OMDB mirror
   (omdb_all_genomes_with_locations.tsv's own genes_aa_file URL for this
   genome -- https://sunagawalab.ethz.ch/.../CARD22-1_SAMN24292811_MAG_00000010.genes.faa.gz,
   5,112 genes total, prodigal headers with real start/end/strand); (b)
   fetched the specific UniProt reference sequences this genome's own
   phaA/phaB/phaE/phaJ/phaZ hits were originally matched to (the
   `best_query` column of PHA_bioprospecting/omdb_search/results/
   pha{A,B,E,J,Z}_unique_targets_with_metadata.tsv, restricted to rows
   where genome == this genome); (c) searched those reference sequences
   against only the ~1,896 genes physically on this genome's 5 phaC-
   bearing scaffolds with pyhmmer's phmmer-style Pipeline.search_seq (a
   fast, local, no-external-binary search -- not re-run automatically by
   this script since it needs the two network fetches above; the
   resulting positions are recorded as constants and can be reproduced
   by anyone re-running steps (a)-(c) against the same URLs/accessions).

3. Structural class of the phaC proteins -- run directly here, fully
   reproducible with no network dependency: the 6 profile HMMs already
   used dataset-wide for the "HMM-supported" evidence tier (Pfam PF07167,
   this project's own 625-column custom HMM, NCBIFam TIGR01838/1839/1836
   for Class I/II/III poly(R)-hydroxyalkanoic acid synthase, PANTHER
   PTHR36837) are searched against these 7 sequences with pyhmmer
   (phac_recovery/hmm/*.hmm, phac_recovery/hmm/extra/*.hmm) -- giving
   real bit scores/e-values per model, not the binary presence-only
   signal available dataset-wide (see section 9.1's correction note on
   why binary tier columns need care). Each TIGRFAM model's own curated
   trusted cutoff (from the HMM file's GA line) is used to say whether a
   target officially passes that class, not just which model scores
   highest.

4. Percent identity of the copies to each other -- the same pairwise
   matrix as (1), shown as a heatmap alongside the tree.

A fifth thing found while cross-checking family assignment, not asked for
directly but load-bearing for how much to trust the two weakest copies:
target_ids ...184878 and ...189118 (the two longest and most divergent
copies, 708aa/736aa, the ones with the weakest/no hits on 4 of 6
structural-class models above) ALSO appear as independent hits in
phaZ_unique_targets_with_metadata.tsv (this genome's own phaZ-family
candidate table) at 44.8%/43.6% identity to their best phaZ reference --
HIGHER than their 39.0%/41.0% identity to their best phaC reference. Both
family assignments are real search hits, not a lookup error; which one is
"correct" is genuinely ambiguous for these two, since phaC and phaZ share
the same alpha/beta-hydrolase fold and catalytic-triad architecture. The
other 5 copies show no such ambiguity (each much more phaC-like than
phaZ-like). PHAC_ZFAMILY_AMBIGUOUS below records this for the plot/TSV.

Usage:
    python figures/scripts/plot_card22_desulfoluna_phac_deep_dive.py

Outputs:
    figures/card22_desulfoluna_phac_deep_dive.png / .pdf
    figures/card22_desulfoluna_phac_paralogs.tsv
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow
import numpy as np
import pyhmmer
from pyhmmer.easel import TextSequence, Alphabet, DigitalSequenceBlock
from pyhmmer.plan7 import HMMFile, Pipeline
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import Bio.Phylo as Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

GENOME = 'CARD22-1_SAMN24292811_MAG_00000010'

# target_id(short) -> (pident to best phaC ref, pident to best phaZ ref) for
# the two copies where family assignment is genuinely ambiguous -- see
# module docstring point 5. From catalytic_domain/phac_sequence_evidence.tsv
# (phaC side) and PHA_bioprospecting/omdb_search/results/
# phaZ_unique_targets_with_metadata.tsv (phaZ side), both already scoped to
# this genome's own target_ids.
PHAC_ZFAMILY_AMBIGUOUS = {
    '184878': (39.0, 44.8),
    '189118': (41.0, 43.6),
}

# ---------------------------------------------------------------------
# 0. one-off, manually-fetched contig data -- see module docstring
#    point 2 for exact provenance. start/end/strand are real prodigal
#    coordinates from this genome's own .genes.faa.gz.
# ---------------------------------------------------------------------
OTHER_PHA_GENES_ON_SCAFFOLDS = [
    # (scaffold, gene_index, start, end, strand, family, best_query_uniprot, note)
    ('scaffold_1', 31, 45892, 46671, -1, 'phaB', 'A0A4Y3PQZ5', ''),
    ('scaffold_1', 56, 74400, 75149, -1, 'phaB', 'A0A0U1P4T1/A0A0W0RRN7', ''),
    ('scaffold_1', 125, 161039, 161800, 1, 'phaB', 'A0A2R8AIW0', ''),
    ('scaffold_1', 236, 302579, 303760, 1, 'phaA', 'A0A0H2ZD32/A0A0S2W1A7', ''),
    ('scaffold_1', 1050, 1303977, 1305851, -1, 'phaA', 'A0A1S6R513', ''),
    ('scaffold_1', 1076, 1337911, 1339116, -1, 'phaA', 'R4KMI7/A0A1G7HQ59/A0A133UVS4/A0A133VJW5', ''),
    ('scaffold_1', 1099, 1371885, 1374122, 1, 'phaA', 'C4XJV8', ''),
    ('scaffold_12', 81, 110267, 111427, 1, 'phaE', 'A0A1W1HIJ3', 'immediately upstream of the phaC copy at idx 82'),
    ('scaffold_12', 83, 112645, 113082, 1, 'phaJ', 'Q142E9/A0ABM7PCC3', 'immediately downstream of the phaC copy at idx 82'),
]
# NOTE: this list is not exhaustive of every phaA/B/E/J/Z gene in the
# genome -- only the ones a direct pyhmmer phmmer-style search of this
# genome's own reference-matched accessions placed ON one of the 5 phaC-
# bearing scaffolds (scaffold_1, 2, 9, 12, 27). scaffold_2, _9, and _27
# had no other-family hits found this way -- their only pha-family genes
# are the phaC copies themselves.

# ---------------------------------------------------------------------
# 1. this genome's 7 verified phaC target_ids (corrected NR100 join)
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
tids = set()
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] == GENOME and row['target_id'] not in bad_targets:
            tids.add(row['target_id'])
tids = sorted(tids)
print(f'{GENOME}: {len(tids)} verified phaC target_ids')

member_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if row[1] in tids:
            member_to_cluster[row[1]] = row[0]

triad = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in tids:
            triad[row['target_id']] = row

ev = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in tids:
            ev[row['target_id']] = row

seqs, positions = {}, {}
cur_id, cur_rep, cur_seq = None, None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in tids:
                seqs[cur_id] = ''.join(cur_seq)
                positions[cur_id] = cur_rep
            parts = line[1:].strip().split(' ', 1)
            cur_id = parts[0]
            cur_rep = parts[1] if len(parts) > 1 else ''
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in tids:
        seqs[cur_id] = ''.join(cur_seq)
        positions[cur_id] = cur_rep


def parse_scaffold_index(rep_str):
    body = rep_str.split('rep=')[-1].split(';')[0]
    if not body.startswith(GENOME + '-'):
        return None, None
    tail = body[len(GENOME) + 1:]
    scaffold, idx = tail.rsplit('_', 1)
    return scaffold, int(idx)


phac_positions = {t: parse_scaffold_index(positions[t]) for t in tids}

# ---------------------------------------------------------------------
# 2. structural class: run the 6 phaC-family HMMs live (pyhmmer, no
#    external binary) -- real bit scores/e-values, not binary presence
# ---------------------------------------------------------------------
alphabet = Alphabet.amino()
digital_seqs = DigitalSequenceBlock(alphabet, [
    TextSequence(name=t.encode(), sequence=seqs[t]).digitize(alphabet) for t in tids
])
HMM_FILES = {
    'PF07167': ROOT / 'phac_recovery/hmm/PF07167.hmm',
    'custom625': ROOT / 'phac_recovery/hmm/phaC_custom.hmm',
    'TIGR01838_ClassI': ROOT / 'phac_recovery/hmm/extra/TIGR01838.hmm',
    'TIGR01839_ClassII': ROOT / 'phac_recovery/hmm/extra/TIGR01839.hmm',
    'TIGR01836_ClassIII': ROOT / 'phac_recovery/hmm/extra/TIGR01836.hmm',
    'PTHR36837': ROOT / 'phac_recovery/hmm/extra/PTHR36837.hmm',
}
CLASS_MODELS = ['TIGR01838_ClassI', 'TIGR01839_ClassII', 'TIGR01836_ClassIII']

hmm_scores = {t: {} for t in tids}
trusted_cutoff = {}
for name, path in HMM_FILES.items():
    with HMMFile(path) as f:
        hmm = next(f)
    trusted_cutoff[name] = hmm.cutoffs.trusted[0] if hmm.cutoffs.trusted else None
    pipeline = Pipeline(alphabet)
    hits = pipeline.search_hmm(hmm, digital_seqs)
    hit_by_name = {}
    for hit in hits:
        hname = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        hit_by_name[hname] = (hit.score, hit.evalue)
    for t in tids:
        hmm_scores[t][name] = hit_by_name.get(t, (None, None))

print('\nStructural class (live pyhmmer search):')
best_class = {}
for t in tids:
    scores = {m: hmm_scores[t][m][0] or 0.0 for m in CLASS_MODELS}
    best = max(scores, key=scores.get)
    passes = trusted_cutoff[best] is not None and scores[best] >= trusted_cutoff[best]
    best_class[t] = (best.replace('TIGR0183_', '').split('_')[-1], passes)
    print(f'  {t[-6:]}: best={best} score={scores[best]:.1f} (trusted cutoff {trusted_cutoff[best]:.1f}) '
          f'{"PASSES" if passes else "below cutoff"}')

# ---------------------------------------------------------------------
# 3. pairwise %identity + NJ tree
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


n = len(tids)
pident_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        p = pident(seqs[tids[i]], seqs[tids[j]])
        pident_matrix[i, j] = pident_matrix[j, i] = p
np.fill_diagonal(pident_matrix, 100.0)

names = [t[-6:] for t in tids]
dm_rows = [[round(1 - pident_matrix[i, j] / 100, 4) for j in range(i + 1)] for i in range(n)]
dm = DistanceMatrix(names=names, matrix=dm_rows)
tree = DistanceTreeConstructor().nj(dm)
tree.ladderize()

# ---------------------------------------------------------------------
# 4. write per-copy TSV
# ---------------------------------------------------------------------
rows_out = []
for t in tids:
    tr = triad.get(t, {})
    e = ev.get(t, {})
    scaf, idx = phac_positions[t]
    cls, passes = best_class[t]
    short = t[-6:]
    zamb = PHAC_ZFAMILY_AMBIGUOUS.get(short)
    rows_out.append({
        'target_id': t, 'protein_length_aa': len(seqs[t]), 'cluster0.7': member_to_cluster.get(t, ''),
        'triad_complete': tr.get('triad_complete', ''), 'pident_to_best_ref': e.get('pident', ''),
        'tcov_to_best_ref': e.get('tcov', ''), 'ref_accession': e.get('ref_accession', ''),
        'best_class_relative': cls, 'passes_class_trusted_cutoff': passes,
        'own_scaffold': scaf or '', 'own_gene_index': idx if idx is not None else '',
        'phaC_vs_phaZ_pident': f'{zamb[0]:.1f} vs {zamb[1]:.1f} (phaZ higher -- ambiguous)' if zamb else '',
    })
out_tsv = OUT / 'card22_desulfoluna_phac_paralogs.tsv'
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)
print(f'\nwrote {out_tsv}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
fig = plt.figure(figsize=(16, 17), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.3, 0.85], hspace=0.5, wspace=0.32)
axTree = fig.add_subplot(gs[0, 0])
axHeat = fig.add_subplot(gs[0, 1])
axContig = fig.add_subplot(gs[1, :])
axClass = fig.add_subplot(gs[2, :])

TEAL, ORANGE, SAGE, MAROON, PURPLE, GOLD, GREY = '#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B', '#7A5FA0', '#C9A227', '#8B8F8C'
FAMILY_COLORS = {'phaA': '#4A7FB5', 'phaB': '#D98E04', 'phaC': SAGE, 'phaE': PURPLE, 'phaJ': MAROON}

# --- Tree panel ---
tid_by_short = {t[-6:]: t for t in tids}


def leaf_color(t):
    return SAGE if triad.get(t, {}).get('triad_complete') == 'True' else ORANGE


label_to_target = {}
for leaf in tree.get_terminals():
    t = tid_by_short[leaf.name]
    cls, passes = best_class[t]
    zflag = ' [phaC/Z ambiguous]' if leaf.name in PHAC_ZFAMILY_AMBIGUOUS else ''
    new_name = f'{leaf.name} ({len(seqs[t])}aa, Cl.{cls}{"*" if passes else ""}){zflag}'
    label_to_target[new_name] = t
    leaf.name = new_name

Phylo.draw(tree, axes=axTree, do_show=False, label_func=lambda c: c.name if c.is_terminal() else '',
           label_colors=lambda lbl: leaf_color(label_to_target.get(lbl, '')))
for text in axTree.texts:
    text.set_fontweight('bold')
    text.set_fontsize(8.7)
axTree.set_title('A. Paralog tree: one clade, or independent classes?', fontsize=12, fontweight='bold', loc='left')
axTree.set_xlabel('substitutions per site (NJ, BLOSUM62 global alignment)', fontsize=8.3)
for spine in axTree.spines.values():
    spine.set_visible(False)
axTree.tick_params(left=False, labelleft=False)
axTree.set_xlim(right=axTree.get_xlim()[1] * 1.55)
sage_p = mpatches.Patch(color=SAGE, label='Triad-complete')
orange_p = mpatches.Patch(color=ORANGE, label='Triad incomplete')
axTree.legend(handles=[sage_p, orange_p], fontsize=7.6, frameon=False, loc='lower right')
axTree.text(0.02, 0.02, '"Cl.X*" = best-fit TIGRFAM class, * = passes official trusted cutoff',
             transform=axTree.transAxes, fontsize=7, color='#5B6E70')

# --- %identity heatmap ---
im = axHeat.imshow(pident_matrix, cmap='YlOrRd', vmin=0, vmax=100)
axHeat.set_xticks(range(n))
axHeat.set_xticklabels(names, fontsize=8, rotation=45, ha='right')
axHeat.set_yticks(range(n))
axHeat.set_yticklabels(names, fontsize=8)
for i in range(n):
    for j in range(n):
        v = pident_matrix[i, j]
        axHeat.text(j, i, f'{v:.0f}', ha='center', va='center', fontsize=7.8,
                     color='white' if v > 55 else '#20302C')
axHeat.set_title('B. Pairwise %identity between all 7 copies', fontsize=12, fontweight='bold', loc='left')
cbar = fig.colorbar(im, ax=axHeat, fraction=0.046, pad=0.04)
cbar.set_label('%identity', fontsize=8.5)
axHeat.text(0.5, -0.22, 'Recent duplicates would cluster >90%; here only two pairs exceed 55% (58%, 67%) --\n'
                       'independent acquisitions/diversification, not one recent duplication.',
             transform=axHeat.transAxes, ha='center', va='top', fontsize=7.6, color='#5B6E70')

# --- Contig map ---
scaffolds_order = ['scaffold_1', 'scaffold_2', 'scaffold_9', 'scaffold_12', 'scaffold_27']
scaffold_genes = {s: [] for s in scaffolds_order}
for t in tids:
    scaf, idx = phac_positions[t]
    if scaf:
        cls, passes = best_class[t]
        zflag = '?' if t[-6:] in PHAC_ZFAMILY_AMBIGUOUS else ''
        scaffold_genes[scaf].append(dict(idx=idx, family='phaC', label=f'phaC{zflag} ({len(seqs[t])}aa, Cl.{cls})',
                                          color=FAMILY_COLORS['phaC'], triad=triad.get(t, {}).get('triad_complete') == 'True'))
for scaf, idx, start, end, strand, fam, acc, note in OTHER_PHA_GENES_ON_SCAFFOLDS:
    scaffold_genes[scaf].append(dict(idx=idx, family=fam, label=fam, color=FAMILY_COLORS[fam], triad=None))

axContig.set_title('C. Contig map: phaC and the other pha-pathway genes on the same scaffolds',
                     fontsize=12, fontweight='bold', loc='left')
axContig.axis('off')
n_scaf = len(scaffolds_order)
slot_h = 1.0 / n_scaf
sub_h = slot_h * 0.5   # leaves a generous gap above/below each track for its
                        # gene-name labels, which render above the arrows and
                        # are not clipped to the inset axes bbox by default
for si, scaf in enumerate(scaffolds_order):
    genes = sorted(scaffold_genes[scaf], key=lambda d: d['idx'])
    sub = axContig.inset_axes([0.03, 1 - (si + 1) * slot_h + slot_h * 0.12, 0.94, sub_h])
    idxs = [g['idx'] for g in genes]
    lo, hi = min(idxs) - 2, max(idxs) + 2
    for g in genes:
        edge = '#20302C' if g['family'] != 'phaC' or g['triad'] else '#B33951'
        lw = 1.6 if g['family'] == 'phaC' else 0.9
        sub.add_patch(FancyArrow(g['idx'] - 0.35, 0, 0.7, 0, width=0.55, head_width=0.8, head_length=0.15,
                                   length_includes_head=True, facecolor=g['color'], edgecolor=edge, linewidth=lw, zorder=3))
        sub.annotate(g['label'], (g['idx'], 0), xytext=(0, 11), textcoords='offset points',
                      ha='center', fontsize=6.6, color='#20302C', rotation=0, clip_on=False)
    sub.set_xlim(lo, hi)
    sub.set_ylim(-1, 1.3)
    sub.axhline(0, color='#C8D2CD', linewidth=0.6, zorder=1)
    sub.set_yticks([])
    sub.set_title(f'{scaf} (gene index)', fontsize=7.8, loc='left', pad=2, color='#3A4A46')
    sub.tick_params(labelsize=7, pad=1)
    for spine in sub.spines.values():
        spine.set_visible(False)

handles = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLORS.items()]
fig.legend(handles=handles, loc='upper center', ncol=5, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, 0.645))

# --- structural-class HMM heatmap ---
model_names = list(HMM_FILES.keys())
score_mat = np.array([[hmm_scores[t][m][0] or 0.0 for m in model_names] for t in tids])
col_max = score_mat.max(axis=0)
col_max[col_max == 0] = 1.0
norm_mat = score_mat / col_max
im2 = axClass.imshow(norm_mat, cmap='YlGnBu', vmin=0, vmax=1, aspect='auto')
axClass.set_xticks(range(len(model_names)))
axClass.set_xticklabels(model_names, fontsize=8.3, rotation=20, ha='right')
axClass.set_yticks(range(n))
axClass.set_yticklabels(names, fontsize=8.3)
for i in range(n):
    for j, m in enumerate(model_names):
        score, evalue = hmm_scores[tids[i]][m]
        cutoff = trusted_cutoff[m]
        if score is None:
            txt = '-'
        else:
            passes = cutoff is not None and score >= cutoff
            txt = f'{score:.0f}{"*" if passes else ""}'
        axClass.text(j, i, txt, ha='center', va='center', fontsize=8,
                       color='white' if norm_mat[i, j] > 0.6 else '#20302C',
                       fontweight='bold' if (score is not None and cutoff and score >= cutoff) else 'normal')
axClass.set_title('D. Structural class: live HMM bit scores per model (color = rank within column; * = passes official trusted cutoff)',
                    fontsize=12, fontweight='bold', loc='left')
fig.subplots_adjust(left=0.07, right=0.97, top=0.94, bottom=0.07)

fig.suptitle('CARD22-1 Desulfoluna MAG (Red Sea coral tissue): 7 phaC copies, one genome',
             fontsize=17, fontweight='bold', y=0.985)
fig.text(0.5, 0.96, f'{GENOME}  ·  99.35% complete, 1.77% contamination  ·  22.18N 39.02E, Goniastrea edwardsi coral tissue  ·  study CARD22-1',
          ha='center', fontsize=9.5, color='#5B6E70')
fig.text(0.5, 0.022, "'phaC?' / '[phaC/Z ambiguous]' = these two copies (184878, 189118) score higher against this genome's own phaZ (depolymerase) "
          'candidates than against phaC references -- family assignment genuinely unresolved, see script docstring point 5.',
          ha='center', fontsize=7.5, color='#5B6E70')

out_path = OUT / 'card22_desulfoluna_phac_deep_dive.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('\nsaved', out_path)
fig.savefig(OUT / 'card22_desulfoluna_phac_deep_dive.pdf', facecolor='white')
print('saved pdf too')
