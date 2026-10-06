"""Phylum x synthase class (I-IV): 100%-stacked bars, one per phylum, of
phaC gene copies (genome x protein pairs) by class. Class type per protein
from catalytic_domain/phac_synthase_class.tsv (best-scoring of the three
NCBIFam class models); class III vs IV split per genome by partner subunit
(phaE -> III, phaR_synthase without phaE -> IV) -- see
classify_phac_synthase_class.py.

Usage:
    python figures/scripts/plot_phac_class_by_phylum.py
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc
import _phylum_colors

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
TOP_N = 10
CLASSES = ['I', 'II', 'III', 'IV', 'III/IV, no partner found', 'unassigned']
# Desaturated greys-to-slate ramp rather than the saturated blue/orange/green/
# amber set this had. In the assembled taxonomy figure the phylum identity is
# carried by the colours from _phylum_colors, in the tree, in the bubble plot and
# in the dot beside each row label here; a second saturated categorical palette
# in the same figure competes with that and makes the panel read as if class and
# phylum were the same kind of variable. Class is an ordered-ish sequence of four
# plus two "unknown" states, so a single-hue ramp suits it and leaves the only
# saturated colour in this panel to the phylum dots.
COLORS = {'I': '#2E4053', 'II': '#5B7793', 'III': '#92ACC4', 'IV': '#C6D4E0',
          'III/IV, no partner found': '#E8ECEF', 'unassigned': '#F5F5F3'}
LIGHT_CLASSES = {'III', 'IV', 'III/IV, no partner found', 'unassigned'}
TEXT_DARK, TEXT_MUTED = '#20302C', '#5B6E70'

sys.path.insert(0, str(ROOT))
from phaatlas.pipeline.pathway_architecture import FAMILY_BAD_QUERIES  # noqa: E402

model_class, confident = {}, {}
with open(ROOT / 'catalytic_domain/phac_synthase_class.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        model_class[row['target_id']] = row['best_model_class']
        confident[row['target_id']] = row['confident'] == 'True'

def uncapped_partner_genomes(family):
    """Genomes carrying a QC-passing hit for a partner family, from the uncapped
    protein-to-genome membership list rather than genome_family_matrix.tsv.

    The matrix is built from a metadata table that lists at most 5 genomes per
    identical protein, so partner evidence is under-recorded exactly the way
    phaC's own count was (PHA_CLEAN_RESULTS.md section 2). The same reference
    exclusion the matrix applies is applied here -- it has to be, and it is the
    dominant filter rather than an afterthought: 94.7% of phaE hits and 97.2% of
    phaR_synthase hits trace to a wrong-gene bait (a Na+/H+ antiporter subunit E
    sharing the Pha name, plus phenylacetate-catabolism proteins; see
    PHA_ALL_FAMILIES_REFERENCE_AUDIT.md). Skipping it would recover 35,733
    phaE genomes that are not phaE.
    """
    path = ROOT / f'{family}_all_genomes_from_nr100_clusters.tsv'
    if not path.exists():
        return set()
    bad_queries = FAMILY_BAD_QUERIES.get(family, frozenset())
    bad_family_targets = set()
    with open(FA / f'{family}_unique_targets_with_metadata.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['best_query'] in bad_queries:
                bad_family_targets.add(row['target_id'])
    genomes = set()
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['target_id'] not in bad_family_targets:
                genomes.add(row['genome'])
    return genomes


phae_extra = uncapped_partner_genomes('phaE')
phar_extra = uncapped_partner_genomes('phaR_synthase')

genome = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if int(row['n_phaC']) > 0:
            name = row['genome']
            genome[name] = (row['gtdb_phylum'],
                            int(row['n_phaE']) > 0 or name in phae_extra,
                            int(row['n_phaR_synthase']) > 0 or name in phar_extra)
print(f'partner evidence: {len(phae_extra):,} phaE and {len(phar_extra):,} phaR_synthase genomes '
      f'from the uncapped lists')

bad = _phac_qc.load_bad_targets()
counts = defaultdict(Counter)
n_conf = n_pairs = 0
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t'); next(r)
    for t, g in r:
        if t in bad or t not in model_class or g not in genome:
            continue
        phylum, has_e, has_r = genome[g]
        mc = model_class[t]
        if mc in ('I', 'II'):
            cls = mc
        elif mc == 'III':
            cls = 'III' if has_e else 'IV' if has_r else 'III/IV, no partner found'
        else:
            cls = 'unassigned'
        counts[phylum][cls] += 1
        n_pairs += 1
        n_conf += confident[t]

total = Counter()
for p, c in counts.items():
    for k, v in c.items():
        total[k] += v
print(f'{n_pairs:,} phaC gene copies in {len(genome):,} genomes; {100 * n_conf / n_pairs:.1f}% above their class model trusted cutoff')
print('overall:', {k: f'{total[k]:,} ({100 * total[k] / n_pairs:.1f}%)' for k in CLASSES})

phyla = [p for p, _ in sorted(counts.items(), key=lambda kv: -sum(kv[1].values())) if p and not p.lower().startswith('unknown')][:TOP_N]
with open(OUT / 'phac_class_by_phylum.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['phylum', 'n_phac_copies'] + CLASSES)
    for p in phyla:
        n = sum(counts[p].values())
        w.writerow([p, n] + [counts[p][c] for c in CLASSES])
        print(f'  {p:20s} n={n:6,}  ' + '  '.join(f'{c.split(",")[0]}:{100 * counts[p][c] / n:.0f}%' for c in CLASSES))

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, ax = plt.subplots(figsize=(8.2, 0.42 * len(phyla) + 1.7), dpi=300)
order = phyla[::-1]
for y, p in enumerate(order):
    n = sum(counts[p].values())
    left = 0
    for c in CLASSES:
        v = 100 * counts[p][c] / n
        if v == 0:
            continue
        ax.barh(y, v, left=left, color=COLORS[c], height=0.68, edgecolor='white', linewidth=1.2)
        if v >= 9:
            ax.text(left + v / 2, y, f'{v:.0f}%', ha='center', va='center', fontsize=8.3,
                    color=TEXT_DARK if c in LIGHT_CLASSES else 'white')
        left += v
    ax.text(101, y, f'n = {n:,}', ha='left', va='center', fontsize=8.3, color=TEXT_MUTED)
ax.set_yticks(range(len(order)))
ax.set_yticklabels(order, fontsize=9.8, color=TEXT_DARK)
# A dot in the phylum's own colour beside each row label, so this panel keys into
# the tree and the bubble plot without repeating their legends.
PHYLUM_DOT_X = -0.016
for y, phylum in enumerate(order):
    ax.scatter([PHYLUM_DOT_X], [y], transform=ax.get_yaxis_transform(), s=66,
               color=_phylum_colors.color_for(phylum), edgecolors='white', linewidths=0.8,
               clip_on=False, zorder=6)
# Tick labels need to clear the dot; without the pad they ran straight through it.
ax.tick_params(axis='y', pad=20)
ax.set_xlim(0, 100)
ax.set_xlabel('Share of phaC gene copies (%)', fontsize=10.5, color=TEXT_DARK)
for side in ('top', 'right', 'left'):
    ax.spines[side].set_visible(False)
ax.spines['bottom'].set_color('#C3C2B7')
ax.tick_params(axis='y', length=0, pad=20)
ax.tick_params(axis='x', colors=TEXT_MUTED, length=3)
handles = [plt.Rectangle((0, 0), 1, 1, facecolor=COLORS[c], edgecolor='none') for c in CLASSES]
labels = ['Class I', 'Class II', 'Class III', 'Class IV', 'III/IV, no partner found', 'unassigned']
fig.legend(handles, labels, loc='upper center', ncol=6, frameon=False, fontsize=8.6, bbox_to_anchor=(0.5, 0.995),
           handlelength=1.1, columnspacing=1.3)
fig.tight_layout(rect=[0, 0, 0.99, 0.94])
fig.savefig(OUT / 'phac_class_by_phylum.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'phac_class_by_phylum.pdf', facecolor='white')
print('saved', OUT / 'phac_class_by_phylum.png')
