"""Substrate-channel profile of each phaC copy in the Desulfoluna multi-copy
MAG, calibrated against solved PHA synthases of known monomer preference.

Reads structure_prediction/phac_pocket_vs_experiment.tsv (run
structure_prediction/anchor_phac_pocket_to_experiment.py first).

Each series is the buried channel next to the catalytic nucleophile, measured
in shells of distance. An extended acyl chain reaches about 1.25 A per carbon
from the thioester sulfur, so the x axis reads as substrate chain length.

The references are the point of the figure. Cupriavidus necator PhaC1 and
Chromobacterium sp. USM2 PhaC are strictly short-chain-length enzymes;
Aeromonas caviae PhaC natively incorporates a C6 monomer. They differ in the
hydrophobicity of the channel at 6-11 A -- exactly where a C6-C8 acyl chain
sits -- far more cleanly than they differ in its volume, so the right-hand
panel is the one that discriminates and the left is context.

Two negative controls are drawn: copy 187467, which has no closed catalytic
triad, and 5T6O, the C319A engineered mutant of the Cupriavidus enzyme. Both
return no channel at any distance, which is what a method measuring the real
active site should do.

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

SOURCE = ROOT / 'structure_prediction/phac_pocket_vs_experiment.tsv'
SHELLS = ['C3-C5', 'C6-C8', 'C9-C12']
SHELL_LABEL = ['0–6 Å\n(C3–C5)', '6–11 Å\n(C6–C8)', '11–16 Å\n(C9–C12)']
TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#B8C0C8', '#EAEDF1'
SCL_COLOR, MCL_COLOR, DEAD = '#8A94A0', '#B3122E', '#C8CDD3'
SHORT_ORGANISM = {'5HZ2': '$\\it{C.\\ necator}$', '5XAV': '$\\it{Chromobacterium}$ sp.'}
# Length, HMM class, and the one model whose confidence is low enough to flag on
# the figure rather than only in the text.
COPY_LABEL = {
    '184878': '184878 · 708 aa · class I',
    '185226': '185226 · 578 aa · class I',
    '185910': '185910 · 577 aa · class III · pLDDT 64',
    '187467': '187467 · 541 aa · no closed triad',
    '187468': '187468 · 380 aa · class III',
    '187880': '187880 · 372 aa · class III',
    '189118': '189118 · 736 aa · class I',
}

style = {c.short: c for c in _deep_dive_copies.copies_for('CARD22-1_SAMN24292811_MAG_00000010')}
rows = list(csv.DictReader(open(SOURCE, newline=''), delimiter='\t'))


def series(row, key):
    return [float(row[f'volume_{s}_A3' if key == 'volume' else f'kd_{s}']) for s in SHELLS]


def spec(row):
    """colour, linestyle, width, marker, zorder, label for one row."""
    if row['kind'] == 'crystal':
        if row['substrate'].startswith('inactive'):
            return DEAD, (0, (2, 3)), 1.4, 'x', 1, f"{row['entry']} — C319A inactive mutant"
        if 'C6' in row['substrate']:
            return MCL_COLOR, (0, (6, 2)), 2.8, 's', 5, f"{row['entry']} — $\\it{{A.\\ caviae}}$, accepts C6"
        return SCL_COLOR, (0, (6, 2)), 2.2, '^', 4, f"{row['entry']} — {SHORT_ORGANISM[row['entry']]}, scl only"
    closed = row['triad_closed'] == 'True'
    colour = style[row['entry']].color if row['entry'] in style else DEAD
    label = COPY_LABEL.get(row['entry'], row['entry'])
    if not closed:
        return DEAD, (0, (2, 3)), 1.4, 'x', 1, COPY_LABEL.get(row['entry'], row['entry'])
    return colour, '-', 2.4, 'o', 3, label


plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0), dpi=300)
x = np.arange(len(SHELLS))

for ax, key, ylabel in ((axes[0], 'volume', 'Connected buried channel volume (Å$^3$)'),
                        (axes[1], 'kd', 'Hydrophobicity of channel lining\n(mean Kyte–Doolittle)')):
    for row in rows:
        colour, ls, lw, marker, z, label = spec(row)
        ax.plot(x, series(row, key), color=colour, linestyle=ls, linewidth=lw, marker=marker,
                markersize=7 if marker != 'x' else 6, markeredgecolor='white',
                markeredgewidth=0.9 if marker != 'x' else 0, zorder=z,
                label=label if ax is axes[0] else None)
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
    ax.set_xlim(-0.3, len(SHELLS) - 0.7)

axes[1].axhline(0, color=RULE, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
axes[1].set_ylim(-0.55, 2.95)
axes[1].annotate('the C6-accepting enzyme stays hydrophobic\nwhere a C6–C8 acyl chain sits',
                 xy=(1.02, 2.42), xytext=(0.62, 2.88), ha='left', va='top', fontsize=8.8,
                 color=MCL_COLOR, linespacing=1.45,
                 arrowprops=dict(arrowstyle='-', color=MCL_COLOR, linewidth=0.9))
axes[1].annotate('strictly short-chain enzymes turn polar here',
                 xy=(1.0, 0.50), xytext=(0.30, -0.42), ha='left', va='bottom', fontsize=8.8,
                 color=TEXT_MUTED,
                 arrowprops=dict(arrowstyle='-', color=SCL_COLOR, linewidth=0.9,
                                 connectionstyle='angle3,angleA=0,angleB=80'))

handles, labels = axes[0].get_legend_handles_labels()
order = sorted(range(len(labels)), key=lambda i: ('—' in labels[i], labels[i]))
fig.legend([handles[i] for i in order], [labels[i] for i in order],
           loc='lower center', ncol=4, frameon=False, fontsize=8.8,
           bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.6)
fig.tight_layout(rect=[0, 0.16, 1, 1])
fig.savefig(ROOT / 'figures/phac_substrate_channel.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/phac_substrate_channel.pdf', facecolor='white')
print('saved', ROOT / 'figures/phac_substrate_channel.png')
