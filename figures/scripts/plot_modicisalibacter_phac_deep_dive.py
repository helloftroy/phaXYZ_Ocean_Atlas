"""Second single-genome deep dive, picked jointly with the user after
comparing two candidates from section 9.3's audit survivors on two
criteria they set: copy count (prefer >5, since the literature maximum
they know of is ~5) and how many sister genomes of the same species exist
for cross-strain comparison.

- `Unknown HK1` (an unnamed lineage, 7 copies in its best genome, 167
  phaC-positive sister genomes across 10 studies) looked promising on
  copy count, but its median copy count across those 167 genomes is only
  2 -- the 7-copy genome is an outlier for the species, not typical, so
  it failed the user's own "median >=3" bar.
- `Modicisalibacter zincidurans` (Halomonadaceae, Gammaproteobacteria --
  the same family as Halomonas bluephagenesis, an already-commercialized
  halotolerant PHA-production chassis) has only 6 sister genomes and none
  with >5 copies, but a direct check of all 6 genomes' 70%-identity
  cluster assignments found THREE phaC paralogs conserved as orthologous
  clusters across literally all 6 independently-sequenced genomes (5
  different studies), with near-identical %identity-to-reference in every
  genome that carries them -- about as clean a "real, conserved species-
  level gene family" signal as this kind of check produces. That result
  is what won the comparison; see PHA_CLEAN_RESULTS.md section 9.5 for
  the full writeup.

Genome: ZHEN20-1_SAMN07748058_MAG_00000118 (western Pacific seawater, near
the Mariana Islands, 11.34N 142.33E) -- the richest of the 6 sister
genomes, carrying a bonus 4th paralog two of the six genomes share.

Four things this script does, extending the CARD22-1 deep-dive script's
method (plot_card22_desulfoluna_phac_deep_dive.py) with one addition:

1. Paralog tree + pairwise %identity for this genome's own 5 copies --
   identical method to the CARD22-1 script (Biopython pairwise BLOSUM62 +
   Bio.Phylo NJ).

2. NEW: cross-genome ortholog conservation -- for the 10 distinct
   target_ids across all 6 sister genomes, shows which 70%-cluster each
   belongs to and which genomes carry it, to visualize the 3 (or 4)
   conserved paralog "slots" directly rather than just reporting the
   count in text.

3. Structural class, live HMM scoring (pyhmmer, same 6 models as the
   CARD22-1 script) -- finds a real disagreement worth reporting: the
   longest copy (738aa, cluster ...753017, confidently phaC by reference
   %identity in every genome that carries it) scores weakly on every
   TIGRFAM class model and has no Pfam PF07167 hit at all, despite that
   strong reference match. Reference-identity confidence and structural-
   HMM confidence do not always agree -- flagged directly rather than
   picking one signal to trust silently.

4. Contig map -- same one-off, documented, manually-fetched-and-recorded
   method as the CARD22-1 script (this genome's own genes.faa.gz from the
   OMDB mirror, plus UniProt sequences for its own phaA/phaB/phaZ hits,
   searched with pyhmmer's phmmer-style Pipeline.search_seq against just
   the genes on the 2 phaC-bearing scaffolds). Finds a genuine, dense
   phaB-phaA-phaB-phaB-phaA-phaC...phaC gene cluster on scaffold_5 --
   the canonical PHA-biosynthesis operon arrangement (phaA+phaB make the
   precursor, phaC polymerizes it) sitting together as one physical unit,
   which is exactly the kind of natural gene cluster that is easiest to
   move into a heterologous host as-is.

One more finding surfaced while checking family assignment, matching what
section 9.4 found in the CARD22-1 genome: two of these 5 copies (both in
the pan-genome-conserved cluster ...783925) score comparably-to-higher
against this genome's own phaZ (depolymerase) candidates than against
phaC references -- PHAC_ZFAMILY_AMBIGUOUS records this, same treatment as
before.

Usage:
    python figures/scripts/plot_modicisalibacter_phac_deep_dive.py

Outputs:
    figures/modicisalibacter_phac_deep_dive.png / .pdf
    figures/modicisalibacter_phac_paralogs.tsv
    figures/modicisalibacter_cross_genome_conservation.tsv
"""
import csv
import sys
from collections import defaultdict
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

GENOME = 'ZHEN20-1_SAMN07748058_MAG_00000118'
SISTER_GENOMES = [
    'DUAR15-1_SAMN06266142_MAG_00000002', 'RSGB23-1_GCF-000731955-V1_GENO_10000001',
    'SANC23-1_SAMEA110646468_MAG_00000065', 'TARA_SAMEA2623562_MAG_00000025',
    'TARA_SAMEA2623564_MAG_00000048', GENOME,
]

# one-off, manually-fetched contig data for GENOME's own phaA/phaB/phaZ
# neighbors -- see module docstring point 4 for exact provenance (same
# method as the CARD22-1 script: this genome's own .genes.faa.gz from the
# OMDB mirror + UniProt reference sequences for its own family hits,
# searched with pyhmmer against just the genes on scaffold_5/scaffold_25).
OTHER_PHA_GENES_ON_SCAFFOLDS = [
    ('scaffold_5', 59, 'phaB', 'multiple phaB refs converge here'),
    ('scaffold_5', 105, 'phaA', ''),
    ('scaffold_5', 108, 'phaB', ''),
    ('scaffold_5', 111, 'phaB', 'multiple phaB refs converge here'),
    ('scaffold_5', 113, 'phaA', 'multiple phaA refs converge here'),
]
# NOTE: not exhaustive of every phaA/B/Z gene genome-wide (n_phaB=12 in
# genome_family_matrix.tsv, only 3 distinct positions found on these 2
# scaffolds) -- only what a direct phmmer search of this genome's own
# reference-matched accessions placed on the 2 phaC-bearing scaffolds.
# scaffold_25 (the genome's third phaC copy) had no other-family hits.

PHAC_ZFAMILY_AMBIGUOUS = {
    '866076': (43.5, 44.5),   # (pident to best phaC ref, pident to best phaZ ref)
    '007752': (43.9, 46.1),
}

# ---------------------------------------------------------------------
# 1. corrected NR100 target_id<->genome join for all 6 sister genomes
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in SISTER_GENOMES and row['target_id'] not in bad_targets:
            genome_targets[row['genome']].add(row['target_id'])

tids = sorted(genome_targets[GENOME])
print(f'{GENOME}: {len(tids)} verified phaC target_ids')
all_sister_targets = set()
for s in genome_targets.values():
    all_sister_targets |= s
print(f'{len(all_sister_targets)} distinct target_ids across all 6 sister genomes')

member_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if row[1] in all_sister_targets:
            member_to_cluster[row[1]] = row[0]

triad = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in all_sister_targets:
            triad[row['target_id']] = row

ev = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in all_sister_targets:
            ev[row['target_id']] = row

seqs, positions = {}, {}
cur_id, cur_rep, cur_seq = None, None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in all_sister_targets:
                seqs[cur_id] = ''.join(cur_seq)
                positions[cur_id] = cur_rep
            parts = line[1:].strip().split(' ', 1)
            cur_id = parts[0]
            cur_rep = parts[1] if len(parts) > 1 else ''
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in all_sister_targets:
        seqs[cur_id] = ''.join(cur_seq)
        positions[cur_id] = cur_rep


def parse_scaffold_index(rep_str):
    body = rep_str.split('rep=')[-1].split(';')[0]
    if '-scaffold_' not in body:
        return None, None
    rep_genome, tail = body.rsplit('-scaffold_', 1)
    if rep_genome != GENOME:
        return None, None
    scaffold, idx = tail.rsplit('_', 1)
    return 'scaffold_' + scaffold, int(idx)


phac_positions = {t: parse_scaffold_index(positions[t]) for t in tids}

# ---------------------------------------------------------------------
# 2. cross-genome ortholog conservation table
# ---------------------------------------------------------------------
cross_rows = []
for g in SISTER_GENOMES:
    for t in sorted(genome_targets[g]):
        e = ev.get(t, {})
        cross_rows.append({
            'genome': g, 'target_id': t, 'cluster0.7': member_to_cluster.get(t, ''),
            'protein_length_aa': len(seqs.get(t, '')), 'pident_to_best_ref': e.get('pident', ''),
            'triad_complete': triad.get(t, {}).get('triad_complete', ''),
        })
cross_out = OUT / 'modicisalibacter_cross_genome_conservation.tsv'
with open(cross_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(cross_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(cross_rows)
print(f'wrote {cross_out}')

clusters_all = sorted(set(r['cluster0.7'] for r in cross_rows), key=lambda c: -sum(1 for r in cross_rows if r['cluster0.7'] == c))
cluster_genome_count = {c: len(set(r['genome'] for r in cross_rows if r['cluster0.7'] == c)) for c in clusters_all}
print('\nOrtholog clusters across the 6 sister genomes:')
for c in clusters_all:
    print(f'  {c[-6:]}: present in {cluster_genome_count[c]}/6 genomes')

# ---------------------------------------------------------------------
# 3. structural class: live 6-HMM search against this genome's 5 copies
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

MEANINGFUL_SCORE_FLOOR = 50  # well below any real class-model signal seen in
# this project's own checked genomes (>100 for a genuinely class-leaning hit,
# often >400 for a confident one) but well above the noise-level 13-33 scores
# a totally unrelated/ambiguous sequence still gets by chance -- below this,
# report "None" rather than naming a class the sequence barely brushes
best_class = {}
for t in tids:
    scores = {m: hmm_scores[t][m][0] or 0.0 for m in CLASS_MODELS}
    best = max(scores, key=scores.get)
    if scores[best] < MEANINGFUL_SCORE_FLOOR:
        best_class[t] = ('None', False)
        continue
    passes = trusted_cutoff[best] is not None and scores[best] >= trusted_cutoff[best]
    best_class[t] = (best.replace('TIGR0183_', '').split('_')[-1], passes)

# ---------------------------------------------------------------------
# 4. pairwise %identity + NJ tree (this genome's own 5 copies)
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
# 5. write per-copy TSV
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
        'n_sister_genomes_sharing_cluster': cluster_genome_count.get(member_to_cluster.get(t, ''), ''),
        'triad_complete': tr.get('triad_complete', ''), 'pident_to_best_ref': e.get('pident', ''),
        'best_class_relative': cls, 'passes_class_trusted_cutoff': passes,
        'own_scaffold': scaf or '', 'own_gene_index': idx if idx is not None else '',
        'phaC_vs_phaZ_pident': f'{zamb[0]:.1f} vs {zamb[1]:.1f} (phaZ comparable/higher -- ambiguous)' if zamb else '',
    })
out_tsv = OUT / 'modicisalibacter_phac_paralogs.tsv'
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)
print(f'wrote {out_tsv}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
fig = plt.figure(figsize=(16, 18), dpi=300)
gs = fig.add_gridspec(4, 2, height_ratios=[1.1, 1.05, 0.85, 0.85], hspace=0.6, wspace=0.32)
axTree = fig.add_subplot(gs[0, 0])
axHeat = fig.add_subplot(gs[0, 1])
axCross = fig.add_subplot(gs[1, :])
axContig = fig.add_subplot(gs[2, :])
axClass = fig.add_subplot(gs[3, :])

TEAL, ORANGE, SAGE, MAROON, PURPLE, GOLD = '#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B', '#7A5FA0', '#C9A227'
FAMILY_COLORS = {'phaA': '#4A7FB5', 'phaB': '#D98E04', 'phaC': SAGE}

# --- A: tree ---
tid_by_short = {t[-6:]: t for t in tids}


def leaf_color(t):
    return SAGE if triad.get(t, {}).get('triad_complete') == 'True' else ORANGE


label_to_target = {}
for leaf in tree.get_terminals():
    t = tid_by_short[leaf.name]
    cls, passes = best_class[t]
    zflag = ' [phaC/Z ambig.]' if leaf.name in PHAC_ZFAMILY_AMBIGUOUS else ''
    nshare = cluster_genome_count.get(member_to_cluster.get(t, ''), '?')
    new_name = f'{leaf.name} ({len(seqs[t])}aa, Cl.{cls}{"*" if passes else ""}, {nshare}/6 genomes){zflag}'
    label_to_target[new_name] = t
    leaf.name = new_name

Phylo.draw(tree, axes=axTree, do_show=False, label_func=lambda c: c.name if c.is_terminal() else '',
           label_colors=lambda lbl: leaf_color(label_to_target.get(lbl, '')))
for text in axTree.texts:
    text.set_fontweight('bold')
    text.set_fontsize(8.2)
axTree.set_title('A. Paralog tree (this genome)', fontsize=12, fontweight='bold', loc='left')
axTree.set_xlabel('substitutions per site (NJ, BLOSUM62)', fontsize=8.3)
for spine in axTree.spines.values():
    spine.set_visible(False)
axTree.tick_params(left=False, labelleft=False)
axTree.set_xlim(right=axTree.get_xlim()[1] * 2.0)
sage_p = mpatches.Patch(color=SAGE, label='Triad-complete')
orange_p = mpatches.Patch(color=ORANGE, label='Triad incomplete')
axTree.legend(handles=[sage_p, orange_p], fontsize=7.6, frameon=False, loc='lower right')

# --- B: %identity heatmap ---
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
axHeat.set_title('B. Pairwise %identity (this genome)', fontsize=12, fontweight='bold', loc='left')
cbar = fig.colorbar(im, ax=axHeat, fraction=0.046, pad=0.04)
cbar.set_label('%identity', fontsize=8.5)

# --- C: cross-genome conservation ---
short_genome_names = {
    'DUAR15-1_SAMN06266142_MAG_00000002': 'DUAR15-1', 'RSGB23-1_GCF-000731955-V1_GENO_10000001': 'RSGB23-1 (isolate)',
    'SANC23-1_SAMEA110646468_MAG_00000065': 'SANC23-1', 'TARA_SAMEA2623562_MAG_00000025': 'TARA-2562',
    'TARA_SAMEA2623564_MAG_00000048': 'TARA-2564', GENOME: 'ZHEN20-1 (this genome)',
}
top_clusters = clusters_all[:4]
mat = np.full((len(SISTER_GENOMES), len(top_clusters)), np.nan)
lenmat = np.full((len(SISTER_GENOMES), len(top_clusters)), np.nan)
for gi, g in enumerate(SISTER_GENOMES):
    for ci, c in enumerate(top_clusters):
        matches = [r for r in cross_rows if r['genome'] == g and r['cluster0.7'] == c]
        if matches:
            mat[gi, ci] = float(matches[0]['pident_to_best_ref']) if matches[0]['pident_to_best_ref'] else np.nan
            lenmat[gi, ci] = matches[0]['protein_length_aa']
masked = np.ma.masked_invalid(mat)
im2 = axCross.imshow(masked, cmap='BuGn', vmin=0, vmax=100, aspect='auto')
axCross.set_xticks(range(len(top_clusters)))
axCross.set_xticklabels([f'paralog {i+1}\n({c[-6:]})' for i, c in enumerate(top_clusters)], fontsize=8.5)
axCross.set_yticks(range(len(SISTER_GENOMES)))
axCross.set_yticklabels([short_genome_names[g] for g in SISTER_GENOMES], fontsize=8.5)
for gi in range(len(SISTER_GENOMES)):
    for ci in range(len(top_clusters)):
        if not np.isnan(mat[gi, ci]):
            axCross.text(ci, gi, f'{mat[gi,ci]:.0f}%\n({lenmat[gi,ci]:.0f}aa)', ha='center', va='center', fontsize=7.6,
                          color='white' if mat[gi, ci] > 55 else '#20302C')
        else:
            axCross.text(ci, gi, '--', ha='center', va='center', fontsize=9, color='#B0B6B2')
axCross.set_title('C. Cross-genome conservation: same paralogs found in all 6 independently-sequenced genomes?',
                    fontsize=12, fontweight='bold', loc='left')
cbar2 = fig.colorbar(im2, ax=axCross, fraction=0.025, pad=0.02)
cbar2.set_label('%identity to reference', fontsize=8)
axCross.text(0.0, -0.22, '5 studies, 1 NCBI isolate reference genome. "--" = this genome lacks that paralog. '
              'Paralogs 1-3 present in all 6 genomes with near-identical %identity/length each time -- a conserved species-level gene family, not noise.',
              transform=axCross.transAxes, fontsize=7.6, color='#5B6E70')

# --- D: contig map ---
scaffolds_order = ['scaffold_5', 'scaffold_25']
scaffold_genes = {s: [] for s in scaffolds_order}
for t in tids:
    scaf, idx = phac_positions[t]
    if scaf:
        cls, passes = best_class[t]
        zflag = '?' if t[-6:] in PHAC_ZFAMILY_AMBIGUOUS else ''
        scaffold_genes[scaf].append(dict(idx=idx, family='phaC', label=f'phaC{zflag} ({len(seqs[t])}aa, Cl.{cls})',
                                          color=FAMILY_COLORS['phaC'], triad=triad.get(t, {}).get('triad_complete') == 'True'))
for scaf, idx, fam, note in OTHER_PHA_GENES_ON_SCAFFOLDS:
    scaffold_genes[scaf].append(dict(idx=idx, family=fam, label=fam, color=FAMILY_COLORS[fam], triad=None))

axContig.set_title('D. Contig map: a real phaB-phaA...phaC gene cluster on scaffold_5',
                     fontsize=12, fontweight='bold', loc='left')
axContig.axis('off')
n_scaf = len(scaffolds_order)
slot_h = 1.0 / n_scaf
sub_h = slot_h * 0.48
for si, scaf in enumerate(scaffolds_order):
    genes = sorted(scaffold_genes[scaf], key=lambda d: d['idx'])
    sub = axContig.inset_axes([0.03, 1 - (si + 1) * slot_h + slot_h * 0.18, 0.94, sub_h])
    idxs = [g['idx'] for g in genes]
    lo, hi = min(idxs) - 2, max(idxs) + 4
    # stagger label height for genes packed close together (label text is
    # much wider than the idx gap between them here) so adjacent labels
    # don't overlap -- close = within 4% of the scaffold's plotted span
    close_thresh = max(1.0, (hi - lo) * 0.10)
    prev_idx = None
    level = 0
    for g in genes:
        if prev_idx is not None and (g['idx'] - prev_idx) < close_thresh:
            level = (level + 1) % 3
        else:
            level = 0
        prev_idx = g['idx']
        y_off = 11 + level * 13
        edge = '#20302C' if g['family'] != 'phaC' or g['triad'] else '#B33951'
        lw = 1.6 if g['family'] == 'phaC' else 0.9
        sub.add_patch(FancyArrow(g['idx'] - 0.35, 0, 0.7, 0, width=0.55, head_width=0.8, head_length=0.15,
                                   length_includes_head=True, facecolor=g['color'], edgecolor=edge, linewidth=lw, zorder=3))
        sub.annotate(g['label'], (g['idx'], 0), xytext=(0, y_off), textcoords='offset points',
                      ha='center', fontsize=6.8, color='#20302C', clip_on=False)
    sub.set_xlim(lo, hi)
    sub.set_ylim(-1, 1.3)
    sub.axhline(0, color='#C8D2CD', linewidth=0.6, zorder=1)
    sub.set_yticks([])
    sub.set_title(f'{scaf} (gene index)', fontsize=7.8, loc='left', pad=2, color='#3A4A46')
    sub.tick_params(labelsize=7, pad=1)
    for spine in sub.spines.values():
        spine.set_visible(False)
handles = [mpatches.Patch(color=c, label=f) for f, c in FAMILY_COLORS.items()]
axContig.legend(handles=handles, loc='upper right', ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(1.0, 1.06))

# --- E: structural-class HMM heatmap ---
model_names = list(HMM_FILES.keys())
score_mat = np.array([[hmm_scores[t][m][0] or 0.0 for m in model_names] for t in tids])
col_max = score_mat.max(axis=0)
col_max[col_max == 0] = 1.0
norm_mat = score_mat / col_max
axClass.imshow(norm_mat, cmap='YlGnBu', vmin=0, vmax=1, aspect='auto')
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
axClass.set_title('E. Structural class: live HMM bit scores (* = passes official trusted cutoff)',
                    fontsize=12, fontweight='bold', loc='left')

fig.suptitle('Modicisalibacter zincidurans: a conserved phaC gene family across 6 independent genomes',
             fontsize=16.5, fontweight='bold', y=0.99)
fig.text(0.5, 0.965, f'{GENOME}  ·  Halomonadaceae (same family as the commercialized PHA-production chassis Halomonas bluephagenesis)  ·  '
          '11.34N 142.33E, western Pacific seawater  ·  study ZHEN20-1',
          ha='center', fontsize=9.3, color='#5B6E70')
fig.text(0.5, 0.022, "'phaC?' / '[phaC/Z ambig.]' = these copies score comparably or higher against this genome's own phaZ "
          '(depolymerase) candidates than against phaC references -- family assignment genuinely unresolved.',
          ha='center', fontsize=7.4, color='#5B6E70')
fig.subplots_adjust(left=0.09, right=0.97, top=0.94, bottom=0.06)

out_path = OUT / 'modicisalibacter_phac_deep_dive.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('\nsaved', out_path)
fig.savefig(OUT / 'modicisalibacter_phac_deep_dive.pdf', facecolor='white')
print('saved pdf too')
