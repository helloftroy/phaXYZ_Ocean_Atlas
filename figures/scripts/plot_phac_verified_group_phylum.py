"""Phylum composition of the verified triad/HMM/pathway groups (see
plot_phac_triad_hmm_pathway_groups.py). Companion to that figure -- tests
whether the taxonomic-bias pattern found pre-audit (Pseudomonadota share
dropping as protein-level evidence drops) survives now that the pathway-
richness tiers are built from verified, audited gene counts.
"""
import pickle
from pathlib import Path
from collections import Counter, defaultdict
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'

genome_group = pickle.load(open('/tmp/phac_verified_triad_hmm_group.pkl', 'rb'))

genome_phylum = {}
with open(ROOT / 'PHA_bioprospecting/omdb_search/results/genome_family_matrix.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        if row['genome'] in genome_group:
            genome_phylum[row['genome']] = row.get('gtdb_phylum') or 'unknown'

ORDER = ['Catalytic triad complete', 'HMM-supported (no triad)', 'No HMM/triad, >=5 other PHA genes',
         'No HMM/triad, 1-4 other PHA genes', 'phaC only']
LABELS = ['Catalytic triad\ncomplete', 'HMM-supported\n(no triad)', '≥5 other\nPHA genes',
          '1-4 other\nPHA genes', 'phaC only']

groups = defaultdict(list)
for g, grp in genome_group.items():
    groups[grp].append(g)

overall = Counter(genome_phylum.values())
TOP_PHYLA = [p for p, _ in overall.most_common(8) if p != 'unknown']
PALETTE = ['#1E6E7A', '#3E8914', '#C2A83E', '#C9622D', '#7A4FA3', '#3E6FC9', '#5B8C7A', '#B4436C', '#9AA098']
colors = dict(zip(TOP_PHYLA + ['Other'], PALETTE))

fracs = {p: [] for p in TOP_PHYLA + ['Other']}
for grp in ORDER:
    gs = groups[grp]
    c = Counter(genome_phylum.get(g, 'unknown') for g in gs)
    n = len(gs)
    for p in TOP_PHYLA + ['Other']:
        # a tier can now be empty (e.g. "phaC only" dropped to 0 genomes after the
        # 2026-09-22 phaC reference-query fix removed exactly the weak, uncorroborated
        # hits that used to populate it) -- 0% rather than a ZeroDivisionError
        val = 100 * (c.get(p, 0) if p in TOP_PHYLA else sum(v for k, v in c.items() if k not in TOP_PHYLA)) / n if n else 0.0
        fracs[p].append(val)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(11, 7.2), dpi=300)

x = range(len(ORDER))
bottom = [0.0] * len(ORDER)
for p in TOP_PHYLA + ['Other']:
    vals = fracs[p]
    ax.bar(x, vals, bottom=bottom, color=colors[p], width=0.62, label=p, zorder=3, edgecolor='white', linewidth=0.4)
    bottom = [b + v for b, v in zip(bottom, vals)]

ax.set_xticks(list(x))
ax.set_xticklabels(LABELS, fontsize=10)
ax.set_ylabel('Share of group genomes (%)', fontsize=12)
ax.set_ylim(0, 100)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=9, frameon=False, title='Phylum (top 8 overall)')

fig.suptitle('Phylum Composition: Verified phaC Evidence Tiers', fontsize=15, fontweight='bold', y=0.99)
fig.text(0.5, 0.925,
          'The pre-audit pattern survives correction: Pseudomonadota share drops as direct protein-level evidence drops.',
          ha='center', fontsize=9.3, color='#5B6E70')

fig.subplots_adjust(left=0.09, right=0.76, top=0.85, bottom=0.11)

out_path = OUT / 'phac_verified_group_phylum_composition.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_verified_group_phylum_composition.pdf', facecolor='white')
