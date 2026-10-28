"""Docking of a chain-length substrate series into each Desulfoluna phaC copy,
shown together with the result that the readout does not pass its own control.

Reads structure_prediction/docking/docking_scores.tsv.

Left: raw Vina affinity against acyl chain length. Every curve falls, because
the ligands gain heavy atoms (13 at C4 to 21 at C12) and score better wherever
they land. What could still be informative is the SLOPE -- a pocket that
accommodates longer chains should fall faster than one that does not.

Right: that slope, as the C4-to-C8 shift, with the span between the two
strictly short-chain reference enzymes shaded. This is the panel that decides
whether any of it can be read. Cupriavidus necator PhaC1 and Chromobacterium
USM2 PhaC are both strictly short-chain-length, yet they sit +0.05 and -1.03
apart; Aeromonas caviae PhaC, which does accept a C6 monomer, reaches -1.35.
Its lead over the nearest same-class reference is 0.32 kcal/mol, smaller than
the 0.72 the two same-class references disagree by. The direction is right and
the separation is not trustworthy.

The figure is drawn this way on purpose. Plotting the per-copy curves without
the reference span would show six confidently separated profiles and invite a
conclusion the data does not support.

Usage:
    python figures/scripts/plot_phac_docking_preference.py
"""
import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _deep_dive_copies

SOURCE = ROOT / 'structure_prediction/docking/docking_scores.tsv'
SCL = ['5HZ2', '5XAV']
MCL = '9KNK'
TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#B8C0C8', '#EAEDF1'
SCL_COLOR, MCL_COLOR, DEAD, BAND = '#8A94A0', '#B3122E', '#C8CDD3', '#E8EAED'

style = {c.short: c for c in _deep_dive_copies.copies_for('CARD22-1_SAMN24292811_MAG_00000010')}
rows = [r for r in csv.DictReader(open(SOURCE, newline=''), delimiter='\t') if r['affinity_kcal_mol']]
raw = defaultdict(dict)
for r in rows:
    raw[r['receptor']][int(r['carbons'])] = float(r['affinity_kcal_mol'])
carbons = sorted({int(r['carbons']) for r in rows})


def spec(name):
    if name in SCL:
        label = {'5HZ2': '$\\it{C.\\ necator}$', '5XAV': '$\\it{Chromobacterium}$ sp.'}[name]
        return SCL_COLOR, (0, (6, 2)), 2.2, '^', 4, f'{name} — {label}, scl only'
    if name == MCL:
        return MCL_COLOR, (0, (6, 2)), 2.8, 's', 5, f'{name} — $\\it{{A.\\ caviae}}$, accepts C6'
    if name == '187467':
        return DEAD, (0, (2, 3)), 1.4, 'x', 1, '187467 — no closed triad'
    return (style[name].color if name in style else DEAD), '-', 2.4, 'o', 3, name


shift = {n: v[8] - v[4] for n, v in raw.items()}
scl_lo, scl_hi = min(shift[n] for n in SCL), max(shift[n] for n in SCL)
noise = st.mean(abs(raw[SCL[0]][c] - raw[SCL[1]][c]) for c in carbons)

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.2), dpi=300,
                         gridspec_kw={'width_ratios': [1.25, 1]})

ax = axes[0]
for name in raw:
    colour, ls, lw, marker, z, label = spec(name)
    ax.plot(carbons, [raw[name][c] for c in carbons], color=colour, linestyle=ls, linewidth=lw,
            marker=marker, markersize=7 if marker != 'x' else 6, markeredgecolor='white',
            markeredgewidth=0.9 if marker != 'x' else 0, zorder=z, label=label)
ax.set_xlabel('Acyl chain length of (R)-3-hydroxyacyl-NAC (carbons)', fontsize=10.2,
              color=TEXT_DARK, labelpad=8)
ax.set_ylabel('AutoDock Vina affinity (kcal/mol)', fontsize=10.2, color=TEXT_DARK, labelpad=8)
ax.set_xticks(carbons)
# Headroom first: the note belongs above the curves, and without extending the
# axis it landed on top of the 5XAV line.
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.55)
ax.text(0.02, 0.985, 'all curves fall with chain length: the C12 ligand has 21 heavy\n'
                     'atoms against C4’s 13, so it scores better wherever it binds',
        transform=ax.transAxes, ha='left', va='top', fontsize=8.7, color=TEXT_MUTED,
        linespacing=1.4)

ax = axes[1]
order = sorted(shift, key=lambda n: shift[n])
ypos = np.arange(len(order))
ax.axvspan(scl_lo, scl_hi, color=BAND, zorder=0)
ax.text((scl_lo + scl_hi) / 2, len(order) + 0.28,
        'span between two enzymes that are\nBOTH strictly short-chain', ha='center', va='top',
        fontsize=8.5, color=TEXT_MUTED, linespacing=1.4)
for y, name in zip(ypos, order):
    colour, ls, lw, marker, z, label = spec(name)
    ax.plot([0, shift[name]], [y, y], color=colour, linewidth=2.0, zorder=2, solid_capstyle='round')
    ax.scatter([shift[name]], [y], s=70, color=colour, marker=marker if marker != 'x' else 'X',
               edgecolors='white', linewidths=0.9, zorder=3)
ax.axvline(0, color=TEXT_DARK, linewidth=1.0, zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels(order, fontsize=9.3, color=TEXT_DARK)
ax.set_ylim(-0.7, len(order) + 0.35)
ax.set_xlabel('C4 → C8 shift in affinity (kcal/mol)\nmore negative = longer chains favoured',
              fontsize=10.2, color=TEXT_DARK, labelpad=8)

for ax in axes:
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(RULE)
    ax.tick_params(colors=TEXT_MUTED, length=3, labelsize=9.3)
axes[1].tick_params(axis='y', length=0)

handles, labels = axes[0].get_legend_handles_labels()
keep = [i for i, l in enumerate(labels) if '—' in l]
fig.legend([handles[i] for i in keep], [labels[i] for i in keep], loc='lower center', ncol=4,
           frameon=False, fontsize=8.8, bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.8)
fig.tight_layout(rect=[0, 0.09, 1, 1])
fig.savefig(ROOT / 'figures/phac_docking_preference.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/phac_docking_preference.pdf', facecolor='white')
print(f'noise floor {noise:.2f}; scl span {scl_lo:+.2f} to {scl_hi:+.2f}; {MCL} {shift[MCL]:+.2f}')
print('saved', ROOT / 'figures/phac_docking_preference.png')
