"""Phylum-controlled habitat effect as a forest plot, built to sit beside
plot_phac_pct_by_ocean_habitat.py as a second column.

The bar chart answers "what fraction of genomes here carry phaC, and how much
of that is explained by which phyla live here". This answers the question that
leaves open: once phylum is controlled for, how big is the habitat effect and
how certain are we? One row per habitat, Cochran-Mantel-Haenszel pooled odds
ratio with a 95% Robins-Breslow-Greenland interval, on a log axis with a line
of no effect at 1.

Pairing rules, so the two figures read as one panel:
  - same rows, same order (ascending raw phaC prevalence, bottom to top), same
    labels, all three taken from _habitat_labels and the stats TSV rather than
    retyped;
  - same two colours, split on the same rule the bar chart uses -- above the
    dataset-wide rate in teal, below it in orange -- so a habitat keeps its
    colour across both panels;
  - the same figure height and margins, so rows line up when placed together.

Input: figures/phac_habitat_phylum_enrichment_test.tsv
  (figures/scripts/phac_habitat_phylum_controlled_test.py -- run it first)

Usage:
    python figures/scripts/plot_phac_habitat_forest.py
"""
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
SOURCE = OUT / 'phac_habitat_phylum_enrichment_test.tsv'

# identical to plot_phac_pct_by_ocean_habitat.py
ACCENT = '#1B7A6E'
LOW_COLOR = '#C9622D'
RULE, GRID = '#B8C0C8', '#E4E8E5'
TEXT_DARK, TEXT_MUTED = '#20302C', '#5B6E70'

if not SOURCE.exists():
    raise SystemExit(f'{SOURCE} is missing -- run figures/scripts/phac_habitat_phylum_controlled_test.py first.')

rows = []
with open(SOURCE, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if not row['cmh_or'] or not row['cmh_or_ci_lo']:
            continue
        rows.append({
            'label': row['habitat_label'], 'n': int(row['n_total']),
            'odds_ratio': float(row['cmh_or']), 'lo': float(row['cmh_or_ci_lo']), 'hi': float(row['cmh_or_ci_hi']),
            'p': float(row['cmh_p']), 'pct_raw': float(row['pct_raw']), 'overall': float(row['overall_pct']),
        })
# ascending raw prevalence, so row i here is row i in the bar chart
rows.sort(key=lambda r: r['pct_raw'])
print(f'{len(rows)} habitats')


def format_p(p):
    """Match the manuscript table: 7.8x10^-85, and a floor below float range."""
    if p == 0:
        return 'p < 10⁻³⁰⁰'
    if p >= 0.001:
        return f'p = {p:.3g}'
    exponent = int(math.floor(math.log10(p)))
    mantissa = p / 10 ** exponent
    superscript = str(exponent).replace('-', '⁻').translate(str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹'))
    return f'p = {mantissa:.1f}×10{superscript}'


# Default is column mode: no habitat names, because the bar chart this sits
# beside already carries them and a second copy reads as a mistake. Row
# geometry matches that figure exactly -- same 0.62 in per row, same 2 in of
# chrome, same tight_layout bottom margin -- so the rows line up when the two
# are placed side by side. --standalone adds the labels back for use on its own.
STANDALONE = '--standalone' in sys.argv
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
fig, ax = plt.subplots(figsize=(5.4 if STANDALONE else 3.6, 0.62 * len(rows) + 2), dpi=300)

y = np.arange(len(rows))
ax.axvline(1.0, color=TEXT_DARK, linewidth=1.0, linestyle=(0, (4, 3)), zorder=2)
for yi, r in enumerate(rows):
    color = ACCENT if r['pct_raw'] >= r['overall'] else LOW_COLOR
    ax.plot([r['lo'], r['hi']], [yi, yi], color=color, linewidth=1.9, solid_capstyle='round', zorder=3)
    for end in (r['lo'], r['hi']):   # end caps, so a short interval is still visible
        ax.plot([end, end], [yi - 0.17, yi + 0.17], color=color, linewidth=1.3, zorder=3)
    ax.scatter([r['odds_ratio']], [yi], s=46, color=color, edgecolors='white', linewidths=1.0, zorder=4)

ax.set_xscale('log')
ax.set_xlim(0.26, 20)
ax.set_xticks([0.5, 1, 2, 5, 10])
ax.set_xticklabels(['0.5', '1', '2', '5', '10'], fontsize=8.5)
ax.set_ylim(-0.75, len(rows) - 0.25)
ax.set_yticks(y)
if STANDALONE:
    ax.set_yticklabels([f"{r['label']}\n(n={r['n']:,})" for r in rows], fontsize=9.3, color=TEXT_DARK)
else:
    ax.set_yticklabels([])
# Short enough to fit the column's own width; the method is named in the caption.
ax.set_xlabel('Odds ratio (95% CI)', fontsize=10, color=TEXT_DARK, labelpad=8)
ax.xaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
for side in ('top', 'right', 'left'):
    ax.spines[side].set_visible(False)
ax.spines['bottom'].set_color(RULE)
ax.tick_params(axis='y', length=0)
ax.tick_params(axis='x', colors=TEXT_MUTED, length=3)

# p-values along the right edge, outside the data, in the row's own colour
for yi, r in enumerate(rows):
    ax.text(1.03, yi, format_p(r['p']), transform=ax.get_yaxis_transform(), ha='left', va='center',
            fontsize=8.2, color=TEXT_MUTED, clip_on=False)

ax.text(1.0, 1.0, 'no effect', transform=ax.get_xaxis_transform(), ha='center', va='bottom',
        fontsize=8.6, color=TEXT_MUTED, style='italic')
ax.text(0.5, -0.062, 'depleted  ←', transform=ax.transAxes, ha='right', va='top',
        fontsize=8.8, color=LOW_COLOR)
ax.text(0.5, -0.062, '  →  enriched', transform=ax.transAxes, ha='left', va='top',
        fontsize=8.8, color=ACCENT)

stem = 'phac_habitat_forest_standalone' if STANDALONE else 'phac_habitat_forest'
fig.tight_layout(rect=[0, 0.045, 0.72, 1])
fig.savefig(OUT / f'{stem}.png', dpi=300, facecolor='white')
fig.savefig(OUT / f'{stem}.pdf', facecolor='white')
print('saved', OUT / f'{stem}.png')
for r in reversed(rows):
    print(f"  {r['label']:34s} OR {r['odds_ratio']:5.2f}  ({r['lo']:.2f}-{r['hi']:.2f})  {format_p(r['p'])}")
