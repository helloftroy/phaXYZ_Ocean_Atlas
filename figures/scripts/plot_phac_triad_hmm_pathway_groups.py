"""Redo of the HMM-support/pathway-context classification for phaC, now
on the fully verified 31,464-genome phaC set (post the 2026-09-22 InterPro
re-verification fix -- see PHA_CLEAN_RESULTS.md; 219 bad reference queries
excluded, up from the earlier 67) and now including the catalytic
Cys-Asp-His triad as the top evidence tier, with pathway context computed
from the VERIFIED other-family counts (genome_family_matrix.tsv, already
regenerated against the corrected exclusion list) rather than the old,
contaminated counts.

Priority order (mutually exclusive, strongest evidence first):
  1. Catalytic triad complete (Cys+Asp+His all present at the expected
     HMM-column positions -- see catalytic_domain/build_triad_table.py)
  2. HMM-supported (>=1 of 6 profile HMMs) but triad not resolvable
  3. Neither HMM nor triad, but >=5 other verified PHA genes present
  4. Neither HMM nor triad, 1-4 other verified PHA genes present
  5. phaC only -- no HMM, no triad, no other PHA gene in the genome

IMPORTANT: genome_targets is filtered against the CURRENT bad-query list
(_phac_qc.load_bad_targets()) every time this runs, rather than trusting
/tmp/verified_phac_genome_targets.pkl to already be filtered -- that pkl
was built against an earlier, smaller exclusion list and going stale
silently is exactly the class of bug that caused rework before. The HMM
hit sets and triad table themselves don't need this treatment: both were
computed against the full 128,199-target pre-audit candidate pool (a
strict superset of any exclusion-list version), so membership-testing
against them stays correct regardless of which targets get filtered out
upstream.
"""
import pickle, csv, sys
from pathlib import Path
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

genome_targets_raw = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
bad_targets = _phac_qc.load_bad_targets()
genome_targets = {}
for g, targets in genome_targets_raw.items():
    kept = targets - bad_targets
    if kept:
        genome_targets[g] = kept
print(f'{len(genome_targets_raw):,} genomes before re-filtering against the current {len(_phac_qc.BAD_QUERIES)}-accession '
      f'exclusion list -> {len(genome_targets):,} genomes after')

hmm_sets = pickle.load(open('/tmp/all_phac_hmm_sets_6models.pkl', 'rb'))
any_hmm = set()
for s in hmm_sets.values():
    any_hmm |= s

triad_targets = set()
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        if row['triad_complete'] == 'True':
            triad_targets.add(row['target_id'])

OTHER_FAMILIES = ['phaA', 'phaB', 'phaE', 'phaR_synthase', 'phaJ', 'phaG', 'phaZ', 'phaY',
                   'phaP', 'phaF', 'phaI', 'phaD', 'phaR_regulator', 'phaQ']

genome_meta = {}
genome_other_count = {}
with open(ROOT / 'PHA_bioprospecting/omdb_search/results/genome_family_matrix.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        g = row['genome']
        if g not in genome_targets:
            continue
        genome_meta[g] = row
        genome_other_count[g] = sum(1 for fam in OTHER_FAMILIES if int(row.get(f'n_{fam}', 0) or 0) > 0)

print(f'{len(genome_meta)} of {len(genome_targets)} verified phaC genomes matched in genome_family_matrix.tsv')

def classify(g):
    targets = genome_targets[g]
    if targets & triad_targets:
        return 'Catalytic triad complete'
    if targets & any_hmm:
        return 'HMM-supported (no triad)'
    n_other = genome_other_count.get(g, 0)
    if n_other >= 5:
        return 'No HMM/triad, >=5 other PHA genes'
    if n_other >= 1:
        return 'No HMM/triad, 1-4 other PHA genes'
    return 'phaC only'

genome_group = {g: classify(g) for g in genome_targets if g in genome_meta}
pickle.dump(genome_group, open('/tmp/phac_verified_triad_hmm_group.pkl', 'wb'))

ORDER = ['Catalytic triad complete', 'HMM-supported (no triad)', 'No HMM/triad, >=5 other PHA genes',
         'No HMM/triad, 1-4 other PHA genes', 'phaC only']
LABELS = ['Catalytic triad\ncomplete', 'HMM-supported\n(no triad)', 'No HMM/triad,\n≥5 other PHA genes',
          'No HMM/triad,\n1-4 other PHA genes', 'phaC only\n(no support)']
COLORS = ['#1E6E7A', '#3E8914', '#7A9B3E', '#C2A83E', '#9E3B3B']

counts = Counter(genome_group.values())
total = sum(counts.values())
values = [counts[k] for k in ORDER]

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(11.5, 7.2), dpi=300)

ax.bar(range(len(ORDER)), values, color=COLORS, width=0.62, zorder=3)
for i, v in enumerate(values):
    ax.text(i, v + total * 0.012, f'{v:,}\n({100*v/total:.1f}%)', ha='center', fontsize=10.5, fontweight='bold', color='#20302C')

ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels(LABELS, fontsize=10)
ax.set_ylabel('Genomes', fontsize=12)
ax.set_ylim(0, max(values) * 1.22)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
ax.set_axisbelow(True)

fig.suptitle('Verified phaC Genomes by Catalytic, HMM, and Pathway Evidence', fontsize=15, fontweight='bold', y=0.975)
fig.text(0.5, 0.925,
          f'All {total:,} genomes with a verified phaC hit (post reference-query audit) -- other-gene counts use the audited, corrected family assignment.',
          ha='center', fontsize=9.3, color='#5B6E70')

footnote = (
    'Priority order (mutually exclusive): catalytic triad > HMM support > pathway richness. Triad = Cys+Asp+His all present at the\n'
    'expected positions via hmmalign-column projection onto phaC_custom.hmm. HMM-supported = hits >=1 of 6 profile HMMs. Other-gene\n'
    'counts exclude phaC itself and use the corrected, per-reference-audited family assignment (see PHA_ALL_FAMILIES_REFERENCE_AUDIT.md).'
)
fig.text(0.5, 0.01, footnote, ha='center', va='bottom', fontsize=7.8, color='#5B6E70')

fig.subplots_adjust(left=0.09, right=0.96, top=0.87, bottom=0.14)

out_path = OUT / 'phac_verified_triad_hmm_pathway_groups.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_verified_triad_hmm_pathway_groups.pdf', facecolor='white')

for k in ORDER:
    print(f'{k}: {counts[k]} ({100*counts[k]/total:.1f}%)')
