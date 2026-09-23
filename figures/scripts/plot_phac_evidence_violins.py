"""Violin-plot view of phaC sequence evidence by evidence tier -- built
because the plain pident-vs-qcov scatter (plot_phac_pident_coverage.py)
overplots into a solid blob at this n and buries the headline finding:
median query coverage stays ~82-84% even in the tiers with no direct
HMM/triad support, so the identity/coverage split by tier is the story,
not a diffuse cloud.
"""
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'

ORDER = ['Catalytic triad complete', 'HMM-supported (no triad)', 'No HMM/triad, >=5 other PHA genes',
         'No HMM/triad, 1-4 other PHA genes', 'phaC only']
LABELS = ['Catalytic triad\ncomplete', 'HMM-supported\n(no triad)', '≥5 other\nPHA genes',
          '1-4 other\nPHA genes', 'phaC only']
COLORS = ['#1E6E7A', '#3E8914', '#7A9B3E', '#C2A83E', '#9E3B3B']

by_tier = defaultdict(lambda: {'pident': [], 'qcov': []})
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        if not row['pident']:
            continue
        t = row['evidence_tier']
        by_tier[t]['pident'].append(float(row['pident']))
        by_tier[t]['qcov'].append(float(row['qcov']) * 100)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 7.5), dpi=300)

def draw_violins(ax, metric, title, ylabel):
    data = [by_tier[t][metric] for t in ORDER]
    # a tier can now be empty (e.g. "phaC only" dropped to 0 genomes after the
    # 2026-09-22 phaC reference-query fix) or too small for a real violin (n<10) --
    # violinplot/boxplot only get the non-empty ones, at their original x positions,
    # so the figure degrades gracefully instead of crashing on an empty array.
    non_empty = [i for i, d in enumerate(data) if d]
    if non_empty:
        parts = ax.violinplot([data[i] for i in non_empty], positions=non_empty, showmedians=False, showextrema=False, widths=0.78)
        for idx, pc in zip(non_empty, parts['bodies']):
            pc.set_facecolor(COLORS[idx])
            pc.set_alpha(0.75)
            pc.set_edgecolor('#2A2A2A')
            pc.set_linewidth(0.6)
        # boxplot overlay for quartiles/median, thin and minimal
        bp = ax.boxplot([data[i] for i in non_empty], positions=non_empty, widths=0.12, patch_artist=True,
                         showfliers=False, medianprops=dict(color='white', linewidth=2),
                         boxprops=dict(facecolor='#2A2A2A', alpha=0.85, edgecolor='none'),
                         whiskerprops=dict(color='#2A2A2A', linewidth=1.3),
                         capprops=dict(color='#2A2A2A', linewidth=1.3))
    for i, d in enumerate(data):
        if not d:
            ax.text(i, 0, 'n=0', va='bottom', ha='center', fontsize=8.5, color='#5B6E70', style='italic')
            continue
        med = np.median(d)
        ax.text(i, med, f' {med:.0f}', va='center', ha='left', fontsize=9, fontweight='bold', color='#111')
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(LABELS, fontsize=9.5)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_ylim(-3, 103)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='#E4E7E2', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)

draw_violins(ax1, 'pident', '% Identity to Best-Matching Audited Reference', '% identity')
draw_violins(ax2, 'qcov', 'Query Coverage to Best-Matching Audited Reference', 'Query coverage (%)')

fig.suptitle('phaC Sequence Evidence by Tier: Identity Drops, Coverage Mostly Doesn\'t', fontsize=15, fontweight='bold', y=0.99)
n_str = ' / '.join(f'{len(by_tier[t]["pident"]):,}' for t in ORDER)
fig.text(0.5, 0.943,
          f'Numbers on each violin are the median. n = {n_str} targets with >=1 hit (left to right).',
          ha='center', fontsize=9.3, color='#5B6E70')

footnote = (
    'Reference set = 1,723 InterPro-confirmed-genuine phaC references (PHA_CLEAN_RESULTS.md §2). Black box = interquartile range, white line = median.\n'
    'The two unresolved tiers are now tiny (n=9 and n=14, versus tens of thousands in the strong-evidence tiers, see PHA_CLEAN_RESULTS.md section 3) --\n'
    'read their violins as a handful of individual remaining cases, not a population estimate. Coverage stays high (~85%) even for these; identity is lower\n'
    'and more variable (median 36-51%) but not below the triad-complete tier\'s own lower spread -- consistent with divergent-but-real matches, not fragments.'
)
fig.text(0.5, 0.01, footnote, ha='center', va='bottom', fontsize=7.9, color='#5B6E70')

fig.subplots_adjust(left=0.07, right=0.97, top=0.86, bottom=0.15, wspace=0.22)

out_path = OUT / 'phac_evidence_violins.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_evidence_violins.pdf', facecolor='white')
