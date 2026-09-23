"""Deep dive on Legionella phaC copy number, following up on section 9.1's
finding that Legionella pneumophila was the single largest species
contributor to the >=5-copy tail (20/239 genomes) -- and that this
population is entirely from RSGB23-1, a bulk NCBI RefSeq isolate-genome
import with no habitat/location metadata (so this script is taxonomy- and
assembly-quality-focused, not ecological).

IMPORTANT CORRECTION to section 9.1's own legitimacy claim, discovered
while building this: that section reported "98.8% of targets are
catalytic-triad-complete" using catalytic_domain/phac_sequence_evidence.tsv's
evidence_tier column. That column is NOT a per-target measurement -- per
catalytic_domain/build_sequence_evidence_table.py, it is each GENOME's
best tier (from phac_verified_triad_hmm_group.pkl) broadcast onto EVERY
target belonging to that genome. For single-copy genomes (the large
majority of the dataset) this is harmless -- one target, one tier, no
difference. For a multi-copy genome it means: if the genome has even ONE
truly triad-complete copy, ALL of its copies get labeled "Catalytic triad
complete" in that table, including copies that individually have zero of
the three catalytic residues. This script instead uses
catalytic_domain/phac_catalytic_triad.tsv's own triad_complete column,
the true per-target measurement, joined via the corrected NR100
target_id<->genome mapping. Section 9.1's own text has been corrected to
match (see PHA_CLEAN_RESULTS.md) -- the true dataset-wide rate is 89.5%
of "tagged complete" targets actually complete (mild inflation, mostly
absorbed by multi-copy genomes), and for the >=5-copy tail specifically,
true per-target completeness is 71.6%, not 98.8%.

Central finding of this script: assembly quality, not biology, mostly
explains this species' copy-number spread. Cross-checked directly against
NCBI's own per-assembly quality metadata (see ASSEMBLY_QC below -- pulled
manually via the NCBI datasets API for a handful of accessions, not
scripted/automated here, since that would need network access this
project's other scripts don't rely on):
  - GCF_000500125.1 (L. pneumophila "W1046", the dataset's single highest
    n_phaC=9 genome) is RefSeq-SUPPRESSED for "many frameshifted
    proteins" -- 440 contigs, N50 16.5kb, 36.9% of its genes are
    pseudogenes, CheckM completeness only 66.75%. Its 9 phaC "copies" are
    mostly short (96-529aa vs. ~590aa full-length), low-target-coverage
    fragments scattered across single-gene contigs -- textbook frameshift
    fragmentation, not real paralogs. Only 2/9 are truly triad-complete.
  - GCF_002082905.2 (L. anisa "FDAARGOS_200", also n_phaC=9) is an FDA-
    ARGOS reference-quality, PacBio+Illumina hybrid assembly: 4 contigs,
    N50 4.2Mb, 99.76% complete, 1.7% pseudogenes, not suppressed. 8/9 of
    its phaC copies are triad-complete and full-length, occupying 9
    DISTINCT 70%-identity clusters (no redundant near-duplicates) --
    genuinely diverse paralogs, not fragmentation.
  - GCF_003205035.1 (L. pneumophila "GC05", n_phaC=6) and GCF_014109805.1
    (Legionella sp. "PC1000", n_phaC=7) are both good/complete assemblies
    (98.9%/99.7% complete, <2.2% pseudogenes, not suppressed) and both
    show 83-86% true per-target triad-completeness -- consistent with the
    two extremes above.
This is a real, checkable pattern, not a guess -- but it was only checked
for 4 of the 27 Legionella genomes with >=5 copies (all from RSGB23-1, an
isolate-genome batch where per-accession NCBI quality metadata actually
exists, unlike this project's MAG-derived genomes). Treat the genus-wide
pathway-architecture correlation below (panel D) with this confound in
mind: an assembly with more frameshift-fragmented phaC copies likely also
has more frameshift-fragmented phaA/phaB/phaZ copies, so part of that
correlation is almost certainly the same artifact acting on every gene
family at once, not coordinated real gene-family expansion.

Usage:
    python figures/scripts/plot_legionella_phac_deep_dive.py

Outputs:
    figures/legionella_phac_deep_dive.png / .pdf
    figures/legionella_phac_copy_matrix.tsv       (all 66 Legionella genomes)
    figures/legionella_phac_top2_paralogs.tsv     (per-copy detail for the two n_phaC=9 genomes)
"""
import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import Bio.Phylo as Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

# manually cross-checked via NCBI datasets API (dataset_report), see module
# docstring -- not re-fetched by this script (no network dependency in this
# project's pipeline scripts otherwise)
ASSEMBLY_QC = {
    'RSGB23-1_GCF-000500125-V1_GENO_10000001': dict(
        accession='GCF_000500125.1', strain='L. pneumophila W1046', n_contigs=440, n50_bp=16550,
        pct_pseudogene=36.9, checkm_complete=66.75, suppressed=True, n_phac_true_triad_frac=2 / 9),
    'RSGB23-1_GCF-002082905-V2_GENO_10000001': dict(
        accession='GCF_002082905.2', strain='L. anisa FDAARGOS_200', n_contigs=4, n50_bp=4234341,
        pct_pseudogene=1.7, checkm_complete=99.76, suppressed=False, n_phac_true_triad_frac=8 / 9),
    'RSGB23-1_GCF-003205035-V1_GENO_10000001': dict(
        accession='GCF_003205035.1', strain='L. pneumophila GC05', n_contigs=112, n50_bp=173439,
        pct_pseudogene=2.1, checkm_complete=98.93, suppressed=False, n_phac_true_triad_frac=5 / 6),
    'RSGB23-1_GCF-014109805-V1_GENO_10000001': dict(
        accession='GCF_014109805.1', strain='Legionella sp. PC1000', n_contigs=4, n50_bp=4081644,
        pct_pseudogene=1.6, checkm_complete=99.7, suppressed=False, n_phac_true_triad_frac=6 / 7),
}
TOP2 = ['RSGB23-1_GCF-000500125-V1_GENO_10000001', 'RSGB23-1_GCF-002082905-V2_GENO_10000001']

# ---------------------------------------------------------------------
# 1. species-by-genome copy-number matrix for the whole genus
# ---------------------------------------------------------------------
genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genomes[row['genome']] = row

leg = {g: r for g, r in genomes.items() if r['gtdb_genus'] == 'Legionella' and int(r['n_phaC']) > 0}
print(f'{len(leg)} Legionella phaC-positive genomes')

fam_cols = ['n_phaA', 'n_phaB', 'n_phaC', 'n_phaJ', 'n_phaZ', 'n_phaY', 'n_phaF', 'n_phaR_regulator']
matrix_out = OUT / 'legionella_phac_copy_matrix.tsv'
with open(matrix_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['genome', 'gtdb_species'] + fam_cols + ['architecture'], delimiter='\t')
    w.writeheader()
    for g, r in sorted(leg.items(), key=lambda kv: -int(kv[1]['n_phaC'])):
        w.writerow({'genome': g, 'gtdb_species': r['gtdb_species'], **{c: r[c] for c in fam_cols}, 'architecture': r['architecture']})
print(f'wrote {matrix_out}')

sp_dist = Counter(r['gtdb_species'] for r in leg.values())
print('species breakdown:', sp_dist.most_common())

# ---------------------------------------------------------------------
# 2. correct per-target join for ALL Legionella (for the true-triad recheck
#    and the architecture correlation's fragmentation caveat)
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
leg_set = set(leg)
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in leg_set and row['target_id'] not in bad_targets:
            genome_targets[row['genome']].add(row['target_id'])

all_leg_targets = set()
for s in genome_targets.values():
    all_leg_targets |= s

triad = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in all_leg_targets:
            triad[row['target_id']] = row

n_true = sum(1 for t in all_leg_targets if triad.get(t, {}).get('triad_complete') == 'True')
print(f'\nLegionella-wide TRUE per-target triad-complete rate: {n_true}/{len(all_leg_targets)} ({100*n_true/len(all_leg_targets):.1f}%)')

# ---------------------------------------------------------------------
# 3. per-copy detail for the two n_phaC=9 genomes: cluster0.7, length,
#    triad, %identity vs reference, and genomic position (scaffold+index
#    from phaC_cluster_sequences.faa's "rep=" field)
# ---------------------------------------------------------------------
member_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        member_to_cluster[row[1]] = row[0]

ev = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ev[row['target_id']] = row

top2_targets = {g: genome_targets[g] for g in TOP2}
all_top2 = set()
for s in top2_targets.values():
    all_top2 |= s

seqs, reps = {}, {}
cur_id, cur_rep, cur_seq = None, None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id and cur_id in all_top2:
                seqs[cur_id] = ''.join(cur_seq)
                reps[cur_id] = cur_rep
            parts = line[1:].strip().split(' ', 1)
            cur_id = parts[0]
            cur_rep = parts[1] if len(parts) > 1 else ''
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id and cur_id in all_top2:
        seqs[cur_id] = ''.join(cur_seq)
        reps[cur_id] = cur_rep


def parse_scaffold_index(rep_str, own_genome):
    # "rep=GENOME-scaffold_N_M;size=S" -- only trust scaffold/index if the
    # rep example is THIS genome itself, not another genome sharing an
    # identical (conserved, cross-strain) sequence
    body = rep_str.split('rep=')[-1].split(';')[0]
    if not body.startswith(own_genome + '-'):
        return None, None
    tail = body[len(own_genome) + 1:]  # "scaffold_N_M"
    parts = tail.rsplit('_', 1)
    scaffold, idx = parts[0], int(parts[1])
    return scaffold, idx


detail_rows = []
for g in TOP2:
    for t in sorted(top2_targets[g]):
        tr = triad.get(t, {})
        e = ev.get(t, {})
        scaffold, idx = parse_scaffold_index(reps.get(t, ''), g)
        detail_rows.append({
            'genome': g, 'target_id': t, 'cluster0.7': member_to_cluster.get(t, ''),
            'protein_length_aa': len(seqs.get(t, '')), 'triad_complete': tr.get('triad_complete', ''),
            'cys_ok': tr.get('cys_ok', ''), 'asp_ok': tr.get('asp_ok', ''), 'his_ok': tr.get('his_ok', ''),
            'pident_to_ref': e.get('pident', ''), 'tcov_to_ref': e.get('tcov', ''),
            'own_scaffold': scaffold or '', 'own_gene_index': idx if idx is not None else '',
        })
detail_out = OUT / 'legionella_phac_top2_paralogs.tsv'
with open(detail_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(detail_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(detail_rows)
print(f'wrote {detail_out}')

n_distinct_clusters = {g: len(set(member_to_cluster.get(t) for t in top2_targets[g])) for g in TOP2}
print(f'\n70%-cluster distinctness: {TOP2[0]} -> {n_distinct_clusters[TOP2[0]]}/9 distinct clusters; '
      f'{TOP2[1]} -> {n_distinct_clusters[TOP2[1]]}/9 distinct clusters')

# ---------------------------------------------------------------------
# 4. pairwise %identity + NJ tree per genome (pure Biopython, no external
#    aligner available in this environment -- global BLOSUM62 alignment)
# ---------------------------------------------------------------------
aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.mode = 'global'


def pident(a, b):
    aln = aligner.align(a, b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    matches = sum(1 for x, y in zip(s1, s2) if x == y and x != '-')
    alen = sum(1 for x, y in zip(s1, s2) if x != '-' and y != '-')
    return 100 * matches / alen if alen else 0.0


trees = {}
for g in TOP2:
    tids = sorted(top2_targets[g])
    names = [t[-6:] for t in tids]
    n = len(tids)
    dm_rows = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                p = pident(seqs[tids[i]], seqs[tids[j]])
                row.append(round(1 - p / 100, 4))
        dm_rows.append(row)
    dm = DistanceMatrix(names=names, matrix=dm_rows)
    tree = DistanceTreeConstructor().nj(dm)
    trees[g] = tree

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.5})
fig = plt.figure(figsize=(16, 15), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[0.85, 1.15, 1.0], hspace=0.55, wspace=0.3)
axA = fig.add_subplot(gs[0, 0])
axQC = fig.add_subplot(gs[0, 1])
axT1 = fig.add_subplot(gs[1, 0])
axT2 = fig.add_subplot(gs[1, 1])
axPos = fig.add_subplot(gs[2, 0])
axArch = fig.add_subplot(gs[2, 1])

TEAL, ORANGE, SAGE, MAROON, GREY = '#1E6E7A', '#C2622D', '#3A6B63', '#9E3B3B', '#8B8F8C'

# Panel A: Legionella copy-count distribution
copy_dist = Counter(int(r['n_phaC']) for r in leg.values())
ks = sorted(copy_dist)
axA.bar([str(k) for k in ks], [copy_dist[k] for k in ks], color=TEAL, zorder=3)
axA.set_xlabel('phaC copies in genome')
axA.set_ylabel('Legionella genomes')
axA.set_title('A. Legionella phaC copy-number distribution', fontsize=11.8, fontweight='bold', loc='left')
axA.text(0.97, 0.95, f'n={len(leg)} phaC+ genomes\n{sp_dist["Legionella pneumophila"]} L. pneumophila\n'
                       f'{len(sp_dist)} species total',
          transform=axA.transAxes, ha='right', va='top', fontsize=8.3, color='#3A4A46',
          bbox=dict(boxstyle='round', facecolor='#F2F5F3', edgecolor='#C8D2CD'))
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.grid(axis='y', color='#E4E8E5', linewidth=0.6, zorder=0)
axA.set_axisbelow(True)

# Panel QC: assembly-quality vs. true per-target triad-completeness (the 4 checked genomes)
qc_genomes = list(ASSEMBLY_QC.keys())
qc_pseudo = [ASSEMBLY_QC[g]['pct_pseudogene'] for g in qc_genomes]
qc_frac = [100 * ASSEMBLY_QC[g]['n_phac_true_triad_frac'] for g in qc_genomes]
qc_labels = [ASSEMBLY_QC[g]['strain'] for g in qc_genomes]
qc_colors = [MAROON if ASSEMBLY_QC[g]['suppressed'] else TEAL for g in qc_genomes]
axQC.scatter(qc_pseudo, qc_frac, s=160, color=qc_colors, edgecolor='#20302C', linewidth=1.0, zorder=3)
label_offsets = {  # hand-placed to avoid the tight cluster of 3 good-quality points colliding
    'L. pneumophila GC05': (8, 14), 'Legionella sp. PC1000': (8, -4), 'L. anisa FDAARGOS_200': (8, -18),
    'L. pneumophila W1046': (-8, 10),
}
for x, y, lbl in zip(qc_pseudo, qc_frac, qc_labels):
    ha = 'right' if label_offsets.get(lbl, (6, 6))[0] < 0 else 'left'
    axQC.annotate(lbl, (x, y), xytext=label_offsets.get(lbl, (6, 6)), textcoords='offset points',
                   fontsize=8, color='#20302C', ha=ha)
axQC.set_xlabel('% of genome\'s genes annotated as pseudogenes (NCBI)')
axQC.set_ylabel('% of its own phaC copies\ntruly triad-complete')
axQC.set_title('QC. Assembly quality predicts copy legitimacy', fontsize=11.8, fontweight='bold', loc='left')
axQC.set_xlim(-2, 40)
axQC.set_ylim(0, 105)
red_patch = mpatches.Patch(color=MAROON, label='RefSeq-suppressed assembly')
teal_patch = mpatches.Patch(color=TEAL, label='Normal assembly')
axQC.legend(handles=[red_patch, teal_patch], fontsize=8, frameon=False, loc='center right')
axQC.spines['top'].set_visible(False)
axQC.spines['right'].set_visible(False)
axQC.grid(color='#E4E8E5', linewidth=0.6, zorder=0)
axQC.set_axisbelow(True)
axQC.text(0.02, 0.03, 'Checked manually via NCBI datasets API for these 4 of 27\nhigh-copy Legionella genomes -- not automated dataset-wide.',
           transform=axQC.transAxes, fontsize=7.3, color='#5B6E70', va='bottom')

# Panels T1/T2: paralog trees for the two n_phaC=9 genomes
tier_by_target = {t: r for t, r in triad.items()}


def status_color(t):
    tr = triad.get(t, {})
    if tr.get('triad_complete') == 'True':
        return SAGE
    return ORANGE


for ax, g, title in [(axT1, TOP2[0], f'T1. {ASSEMBLY_QC[TOP2[0]]["strain"]} (RefSeq-suppressed)'),
                      (axT2, TOP2[1], f'T2. {ASSEMBLY_QC[TOP2[1]]["strain"]} (reference-quality)')]:
    tree = trees[g]
    tree.ladderize()
    tid_by_short = {t[-6:]: t for t in top2_targets[g]}
    label_to_target = {}
    for leaf in tree.get_terminals():
        t = tid_by_short[leaf.name]
        new_name = f'{leaf.name} ({len(seqs[t])}aa)'
        label_to_target[new_name] = t
        leaf.name = new_name

    def color_for_label(label, mapping=label_to_target):
        t = mapping.get(label)
        return status_color(t) if t else 'black'

    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: c.name if c.is_terminal() else '',
               label_colors=color_for_label)
    for text in ax.texts:
        text.set_fontweight('bold')
        text.set_fontsize(9)
    ax.set_title(title, fontsize=11.5, fontweight='bold', loc='left')
    ax.set_xlabel('substitutions per site (NJ, BLOSUM62 global alignment)', fontsize=8.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.set_xlim(right=ax.get_xlim()[1] * 1.35)  # room for tip labels before the panel edge

sage_patch = mpatches.Patch(color=SAGE, label='Triad-complete')
orange_patch = mpatches.Patch(color=ORANGE, label='Triad incomplete/fragment')
fig.legend(handles=[sage_patch, orange_patch], fontsize=9, frameon=False, ncol=2,
           loc='upper center', bbox_to_anchor=(0.5, 0.635))

# Panel Pos: genomic position (scaffold, gene index) for both genomes --
# two stacked strip plots (own scaffold vs. gene index), one per genome
axPos.axis('off')
sub1 = axPos.inset_axes([0.0, 0.55, 1.0, 0.34])
sub2 = axPos.inset_axes([0.0, 0.02, 1.0, 0.34])
for ax_s, g, color, title in [(sub1, TOP2[0], MAROON, ASSEMBLY_QC[TOP2[0]]['strain']),
                                (sub2, TOP2[1], TEAL, ASSEMBLY_QC[TOP2[1]]['strain'])]:
    rows_g = [r for r in detail_rows if r['genome'] == g and r['own_scaffold']]
    scafs = sorted(set(r['own_scaffold'] for r in rows_g))
    scaf_y = {s: i for i, s in enumerate(scafs)}
    for r in rows_g:
        y = scaf_y[r['own_scaffold']]
        c = SAGE if r['triad_complete'] == 'True' else ORANGE
        ax_s.scatter(r['own_gene_index'], y, s=90, color=c, edgecolor='#20302C', linewidth=0.8, zorder=3)
        ax_s.annotate(f"{r['protein_length_aa']}aa", (r['own_gene_index'], y), xytext=(0, 7),
                       textcoords='offset points', ha='center', fontsize=6.8, color='#3A4A46')
    ax_s.set_yticks(range(len(scafs)))
    ax_s.set_yticklabels([s.replace('scaffold_', 'scf.') for s in scafs], fontsize=7.5)
    ax_s.set_ylim(-0.7, len(scafs) - 0.3)
    ax_s.set_xlabel('gene index on its own scaffold', fontsize=8)
    ax_s.set_title(title, fontsize=9.5, fontweight='bold', loc='left', color=color, pad=14)
    ax_s.spines['top'].set_visible(False)
    ax_s.spines['right'].set_visible(False)
    ax_s.grid(axis='x', color='#E4E8E5', linewidth=0.5, zorder=0)
    ax_s.set_axisbelow(True)
axPos.set_title('Pos. Genomic position: same scaffold or scattered?', fontsize=11.8, fontweight='bold', loc='left', y=1.18)

# Panel Arch: pathway architecture, high- vs low-copy Legionella
high = [r for r in leg.values() if int(r['n_phaC']) >= 5]
low = [r for r in leg.values() if int(r['n_phaC']) <= 2]
arch_cols = ['n_phaA', 'n_phaB', 'n_phaJ', 'n_phaZ', 'n_phaY', 'n_phaF', 'n_phaR_regulator']
high_means = [np.mean([int(r[c]) for r in high]) for c in arch_cols]
low_means = [np.mean([int(r[c]) for r in low]) for c in arch_cols]
y = np.arange(len(arch_cols))
h = 0.36
axArch.barh(y + h / 2, high_means, height=h, color=MAROON, label=f'n_phaC>=5 (n={len(high)})', zorder=3)
axArch.barh(y - h / 2, low_means, height=h, color=TEAL, label=f'n_phaC<=2 (n={len(low)})', zorder=3)
axArch.set_yticks(y)
axArch.set_yticklabels([c.replace('n_pha', 'pha') for c in arch_cols])
axArch.set_xlabel('Mean gene count per genome')
axArch.set_title('Arch. Other PHA genes scale with phaC copy number too', fontsize=11.8, fontweight='bold', loc='left')
axArch.legend(fontsize=8.3, frameon=False)
axArch.spines['top'].set_visible(False)
axArch.spines['right'].set_visible(False)
axArch.spines['left'].set_visible(False)
axArch.tick_params(axis='y', length=0)
axArch.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axArch.set_axisbelow(True)
axArch.text(0.98, 0.03, 'Caveat: likely partly the same fragmentation\nartifact acting on every gene family at once (see QC panel).',
             transform=axArch.transAxes, ha='right', fontsize=7.3, color='#5B6E70')

fig.suptitle('Legionella phaC: a genus-wide signal mostly explained by assembly quality',
             fontsize=17, fontweight='bold', y=0.99)
fig.text(0.5, 0.965,
          f'96% of L. pneumophila genomes (49/51) carry phaC, most with multiple copies -- but per-copy '
          f'legitimacy tracks assembly quality (QC panel), not a uniform biological signal.',
          ha='center', fontsize=9.7, color='#5B6E70')
fig.subplots_adjust(left=0.08, right=0.97, top=0.925, bottom=0.055)

out_path = OUT / 'legionella_phac_deep_dive.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('\nsaved', out_path)
fig.savefig(OUT / 'legionella_phac_deep_dive.pdf', facecolor='white')
print('saved pdf too')
