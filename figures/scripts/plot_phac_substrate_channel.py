"""Substrate-channel profile of each phaC copy in the Desulfoluna multi-copy MAG.

Reads structure_prediction/phac_substrate_pocket.tsv (run
structure_prediction/analyze_phac_substrate_pocket.py first) and draws, per
copy, the buried channel volume and its hydrophobicity at increasing distance
from the catalytic nucleophile. Shell boundaries correspond to the reach of an
extended acyl chain at roughly 1.25 A per carbon, so the x axis reads as
substrate chain length.

187467 is kept in as a negative control, not padding: it is the copy with no
closed catalytic triad, and it returns no connected channel at any distance.
A method that reports a plausible pocket for it would not be measuring the
active site.

Usage:
    python figures/scripts/plot_phac_substrate_channel.py
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _deep_dive_copies

SOURCE = ROOT / 'structure_prediction/phac_substrate_pocket.tsv'
SHELLS = ['C3-C5', 'C6-C8', 'C9-C12']
SHELL_LABEL = ['0–6 Å\n(C3–C5)', '6–11 Å\n(C6–C8)', '11–16 Å\n(C9–C12)']
TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#B8C0C8', '#EAEDF1'
DEAD = '#9AA5B1'

style = {c.short: c for c in _deep_dive_copies.copies_for('CARD22-1_SAMN24292811_MAG_00000010')}
rows = list(csv.DictReader(open(SOURCE, newline=''), delimiter='\t'))
print(f'{len(rows)} copies')

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), dpi=300)
x = np.arange(len(SHELLS))

for ax, key, ylabel in ((axes[0], 'volume', 'Connected buried channel volume (Å$^3$)'),
                        (axes[1], 'kd', 'Mean Kyte–Doolittle hydropathy of lining residues')):
    for row in rows:
        closed = row['triad_closed'] == 'True'
        colour = style[row['copy']].color if row['copy'] in style else DEAD
        values = [float(row[f'{key}_{s}_A3' if key == 'volume' else f'{key}_{s}']) for s in SHELLS]
        ax.plot(x, values, color=colour, linewidth=2.4 if closed else 1.6,
                linestyle='-' if closed else (0, (4, 3)),
                marker='o' if closed else 'x', markersize=7, markeredgecolor='white',
                markeredgewidth=0.9 if closed else 0, zorder=3 if closed else 2,
                label=f"{row['copy']} · {row['length_aa']} aa · class {row['model_class']}"
                      + ('' if closed else '  (no closed triad)'))
    ax.set_xticks(x)
    ax.set_xticklabels(SHELL_LABEL, fontsize=9.3, color=TEXT_DARK)
    ax.set_xlabel('Distance from catalytic nucleophile', fontsize=10.5, color=TEXT_DARK, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=10.3, color=TEXT_DARK, labelpad=8)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(RULE)
    ax.tick_params(colors=TEXT_MUTED, length=3, labelsize=9.3)
    ax.set_xlim(-0.25, len(SHELLS) - 0.75)

axes[1].axhline(0, color=RULE, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
axes[1].text(-0.18, 0.04, 'hydrophobic above', ha='left', va='bottom',
             fontsize=8.4, color=TEXT_MUTED)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, fontsize=9.0,
           bbox_to_anchor=(0.5, 0.0), handlelength=2.0, columnspacing=2.0)
fig.tight_layout(rect=[0, 0.14, 1, 1])
fig.savefig(ROOT / 'figures/phac_substrate_channel.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/phac_substrate_channel.pdf', facecolor='white')
print('saved', ROOT / 'figures/phac_substrate_channel.png')
