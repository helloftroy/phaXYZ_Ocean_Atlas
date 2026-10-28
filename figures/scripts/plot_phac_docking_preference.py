"""Chain-length preference of each Desulfoluna phaC copy from docking.

Reads structure_prediction/docking/docking_preference.tsv (run
structure_prediction/docking/analyse_docking.py first).

Left panel is the raw Vina affinity and exists to show the problem rather than
the result: every protein's curve falls with chain length because the ligands
get bigger, not because the pockets prefer them. Right panel subtracts, for
each ligand separately, the mean affinity the two strictly short-chain-length
reference enzymes achieve with that same ligand. Size cancels within a ligand,
so what is left is preference relative to enzymes whose preference is known.

Below zero on the right means the protein accommodates that chain length
better than Cupriavidus necator and Chromobacterium USM2 PhaC do.

Usage:
    python figures/scripts/plot_phac_docking_preference.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _deep_dive_copies

SOURCE = ROOT / 'structure_prediction/docking/docking_preference.tsv'
TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#B8C0C8', '#EAEDF1'
SCL_COLOR, MCL_COLOR, DEAD = '#8A94A0', '#B3122E', '#C8CDD3'

style = {c.short: c for c in _deep_dive_copies.copies_for('CARD22-1_SAMN24292811_MAG_00000010')}
rows = list(csv.DictReader(open(SOURCE, newline=''), delimiter='\t'))

series = defaultdict(dict)
meta = {}
for r in rows:
    series[r['receptor']][int(r['carbons'])] = r
    meta[r['receptor']] = r


def spec(name):
    r = meta[name]
    if r['kind'] == 'crystal':
        if 'C6' in r['substrate']:
            return MCL_COLOR, (0, (6, 2)), 2.8, 's', 5, f"{name} — $\\it{{A.\\ caviae}}$, accepts C6"
        short = {'5HZ2': '$\\it{C.\\ necator}$', '5XAV': '$\\it{Chromobacterium}$ sp.'}[name]
        return SCL_COLOR, (0, (6, 2)), 2.2, '^', 4, f'{name} — {short}, scl only'
    if name == '187467':
        return DEAD, (0, (2, 3)), 1.4, 'x', 1, '187467 — no closed triad'
    colour = style[name].color if name in style else DEAD
    return colour, '-', 2.4, 'o', 3, r['label']


plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0), dpi=300)

panels = [(axes[0], 'affinity_kcal_mol', 'Raw Vina affinity (kcal/mol)'),
          (axes[1], 'vs_scl_reference', 'Affinity relative to the short-chain\nreference enzymes (kcal/mol)')]
for ax, key, ylabel in panels:
    for name in series:
        colour, ls, lw, marker, z, label = spec(name)
        carbons = sorted(series[name])
        values = [float(series[name][c][key]) for c in carbons if series[name][c][key] != '']
        xs = [c for c in carbons if series[name][c][key] != '']
        ax.plot(xs, values, color=colour, linestyle=ls, linewidth=lw, marker=marker,
                markersize=7 if marker != 'x' else 6, markeredgecolor='white',
                markeredgewidth=0.9 if marker != 'x' else 0, zorder=z,
                label=label if ax is axes[0] else None)
    ax.set_xlabel('Acyl chain length of (R)-3-hydroxyacyl-NAC substrate (carbons)',
                  fontsize=10.2, color=TEXT_DARK, labelpad=8)
    ax.set_ylabel(ylabel, fontsize=10.2, color=TEXT_DARK, labelpad=8)
    ax.set_xticks(sorted({int(r['carbons']) for r in rows}))
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(RULE)
    ax.tick_params(colors=TEXT_MUTED, length=3, labelsize=9.3)

axes[0].text(0.5, 0.02, 'every curve falls with chain length:\nbigger ligands score better regardless of fit',
             transform=axes[0].transAxes, ha='center', va='bottom', fontsize=8.8,
             color=TEXT_MUTED, linespacing=1.4)
axes[1].axhline(0, color=TEXT_DARK, linewidth=1.1, linestyle=(0, (4, 3)), zorder=2)
axes[1].text(axes[1].get_xlim()[1], 0.06, 'short-chain references', ha='right', va='bottom',
             fontsize=8.6, color=TEXT_MUTED, transform=axes[1].get_yaxis_transform(which='grid'))

handles, labels = axes[0].get_legend_handles_labels()
order = sorted(range(len(labels)), key=lambda i: ('—' in labels[i], labels[i]))
fig.legend([handles[i] for i in order], [labels[i] for i in order], loc='lower center', ncol=4,
           frameon=False, fontsize=8.8, bbox_to_anchor=(0.5, 0.0), handlelength=2.2, columnspacing=1.6)
fig.tight_layout(rect=[0, 0.16, 1, 1])
fig.savefig(ROOT / 'figures/phac_docking_preference.png', dpi=300, facecolor='white')
fig.savefig(ROOT / 'figures/phac_docking_preference.pdf', facecolor='white')
print('saved', ROOT / 'figures/phac_docking_preference.png')
