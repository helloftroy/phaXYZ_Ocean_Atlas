"""Manuscript version of the sequence-vs-structure figure (section 9.8).

Same data and same claim as plot_seq_vs_structure_hero.py -- one point per
all_phac_dedup representative, sequence identity to its best-matching
reference against structural similarity (qTM) to that same reference, split
by catalytic-triad completeness -- redrawn for print. That script's own
output is read back from figures/seq_vs_structure_hero.tsv rather than
rebuilding the join, so the two figures cannot drift apart; run it first if
the underlying tables change.

What changed and why, since the data is identical:

  the bubbles are gone.  Size encoded how many genomes carry each
    representative's paralog cluster. With 15,867 points in one panel it
    cost far more than it bought: large bubbles swallowed their neighbours,
    the size legend needed a second boxed legend inside the plot, and the
    cluster-size question is answered properly by its own figure
    (phac_cluster_occupancy_spectrum). Points are now uniform and small.
  the cloud is readable.  Points are drawn small, translucent, rasterized
    (so the PDF stays light) and in random order. Random order matters: any
    fixed draw order puts one category wholly on top of the other and makes
    whichever is drawn last look more abundant than it is.
  the signal is drawn, not inferred.  A running median with an
    interquartile ribbon per group now runs across the cloud. The trend was
    previously left for the reader to find by eye in an overplotted mass.
  a marginal panel carries the distribution.  qTM densities per group, on
    the shared y axis, which is where the separation between the two groups
    actually lives.
  nothing is written inside the panel.  No title, no subtitle, no boxed
    legends over data. The caption belongs in the manuscript.

Usage:
    python figures/scripts/plot_seq_vs_structure_publication.py

Outputs:
    figures/seq_vs_structure.png / .pdf
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
SOURCE = OUT / 'seq_vs_structure_hero.tsv'

# Blue from the sequential ramp the other manuscript figures use, against the
# warm terracotta this comparison has always used. Blue/orange is the pair that
# survives every common form of colour-vision deficiency and greyscale print,
# which matters more here than hue continuity with the single-series figures.
TRIAD_COLOR = '#1c5cab'
NO_TRIAD_COLOR = '#C2622D'
BAR_LIGHT, BAR_DARK = '#9ec5f4', '#1c5cab'
GRID, TWILIGHT, RULE = '#EDEFF2', '#F4F6F8', '#B8C0C8'
TEXT_DARK, TEXT_MUTED = '#1E2630', '#5B6670'

BINS = [(0, 30, '<30%'), (30, 40, '30–40%'), (40, 50, '40–50%'), (50, 70, '50–70%'), (70, 100.0001, '≥70%')]
TWILIGHT_MAX = 30

if not SOURCE.exists():
    raise SystemExit(f'{SOURCE} is missing -- run figures/scripts/plot_seq_vs_structure_hero.py first.')

rows = []
with open(SOURCE, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        rows.append({'pident': float(row['pident']), 'qtm': 100 * float(row['qtmscore']),
                     'triad': row['triad_complete'] == 'True'})
print(f'{len(rows):,} representatives')

x = np.array([r['pident'] for r in rows])
y = np.array([r['qtm'] for r in rows])
triad = np.array([r['triad'] for r in rows])
n_triad, n_no_triad = int(triad.sum()), int((~triad).sum())

pct_real = 100 * np.mean(y >= 50)
twilight = x < TWILIGHT_MAX
print(f'{pct_real:.1f}% at qTM>=50; twilight zone n={twilight.sum():,}, '
      f'{100 * np.mean(y[twilight] >= 85):.1f}% at qTM>=85, {100 * np.mean(y[twilight] >= 50):.1f}% at qTM>=50')

bin_stats = []
for lo, hi, label in BINS:
    sub = y[(x >= lo) & (x < hi)]
    bin_stats.append({'label': label, 'n': len(sub),
                      'p50': 100 * np.mean(sub >= 50), 'p85': 100 * np.mean(sub >= 85)})
    print(f"  {label:7s} n={len(sub):6,d}  qTM>=50: {bin_stats[-1]['p50']:5.1f}%  qTM>=85: {bin_stats[-1]['p85']:5.1f}%")

# ---------------------------------------------------------------- figure
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
fig = plt.figure(figsize=(9.6, 8.4), dpi=300)
gs = fig.add_gridspec(2, 2, width_ratios=[1, 0.19], height_ratios=[0.30, 1],
                      hspace=0.17, wspace=0.035, left=0.085, right=0.975, top=0.945, bottom=0.115)
ax_bin = fig.add_subplot(gs[0, 0])
ax = fig.add_subplot(gs[1, 0])
ax_marg = fig.add_subplot(gs[1, 1], sharey=ax)
fig.add_subplot(gs[0, 1]).axis('off')

# ---- top: structural support by sequence-identity bin
positions = np.arange(len(bin_stats))
ax_bin.grid(True, axis='y', color=GRID, linewidth=0.8, zorder=0)
ax_bin.set_axisbelow(True)
for pos, b in zip(positions, bin_stats):
    ax_bin.bar(pos, b['p50'], width=0.56, color=BAR_LIGHT, zorder=2, edgecolor='none')
    ax_bin.bar(pos, b['p85'], width=0.28, color=BAR_DARK, zorder=3, edgecolor='none')
    ax_bin.text(pos, b['p50'] + 4, f"{b['p50']:.0f}%", ha='center', va='bottom', fontsize=8.6, color=TEXT_MUTED)
    # A short dark bar can hold no label at all: white on it is fine only while
    # the bar is tall enough to contain the text, and dark text on dark navy is
    # unreadable. Below that height the label goes just above the bar, where the
    # pale bar behind it gives the contrast instead.
    if b['p85'] > 25:
        ax_bin.text(pos, b['p85'] - 16, f"{b['p85']:.0f}%", ha='center', va='bottom', fontsize=8.6,
                    color='white', fontweight='bold')
    else:
        ax_bin.text(pos, b['p85'] + 3, f"{b['p85']:.0f}%", ha='center', va='bottom', fontsize=8.6,
                    color=BAR_DARK, fontweight='bold')
    ax_bin.text(pos, -9, f"{b['label']}\nn={b['n']:,}", ha='center', va='top', fontsize=8.4, color=TEXT_MUTED,
                linespacing=1.4)
ax_bin.set_xlim(-0.62, len(bin_stats) - 0.38)
ax_bin.set_ylim(0, 108)
ax_bin.set_xticks([])
ax_bin.set_yticks([0, 50, 100])
ax_bin.set_yticklabels(['0', '50', '100'], fontsize=8.6)
ax_bin.set_ylabel('% of bin', fontsize=9, color=TEXT_DARK)
for spine in ('top', 'right', 'bottom'):
    ax_bin.spines[spine].set_visible(False)
ax_bin.spines['left'].set_color(RULE)
ax_bin.tick_params(axis='y', length=0, colors=TEXT_MUTED)
ax_bin.legend(handles=[mpatches.Patch(color=BAR_LIGHT, label='qTM ≥ 50 (structurally real)'),
                       mpatches.Patch(color=BAR_DARK, label='qTM ≥ 85 (confident)')],
              loc='lower right', bbox_to_anchor=(1.0, 1.0), fontsize=8.6, frameon=False, ncol=2,
              handlelength=1.3, columnspacing=1.4)

# ---- main panel
ax.axvspan(18, TWILIGHT_MAX, color=TWILIGHT, zorder=0, linewidth=0)
ax.axvline(TWILIGHT_MAX, color=RULE, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
for threshold, label in ((85, 'qTM ≥ 85'), (50, 'qTM ≥ 50')):
    ax.axhline(threshold, color=RULE, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
    ax.text(99.4, threshold + 1.0, label, ha='right', va='bottom', fontsize=8.4, color=TEXT_MUTED,
            bbox=dict(boxstyle='round,pad=0.22', facecolor='white', edgecolor='none', alpha=0.85))

# Random draw order, so neither category is systematically painted over the
# other; a fixed order makes whichever is drawn last look more abundant.
order = np.random.default_rng(42).permutation(len(x))
ax.scatter(x[order], y[order], s=3.2, linewidth=0, alpha=0.22, rasterized=True, zorder=2,
           c=np.where(triad[order], TRIAD_COLOR, NO_TRIAD_COLOR))


def running_summary(mask, width=6, step=2):
    """Median and IQR of qTM in a sliding window of sequence identity. Windows
    with fewer than 25 points are dropped rather than plotted as a noisy
    median -- the low-identity end is sparse and would otherwise end in a spike
    that looks like signal."""
    centres, med, lo, hi = [], [], [], []
    for centre in np.arange(width / 2, 100 - width / 2 + step, step):
        window = mask & (x >= centre - width / 2) & (x < centre + width / 2)
        if window.sum() < 25:
            continue
        values = y[window]
        centres.append(centre)
        med.append(np.median(values))
        lo.append(np.percentile(values, 25))
        hi.append(np.percentile(values, 75))
    return np.array(centres), np.array(med), np.array(lo), np.array(hi)


for mask, color, label in ((triad, TRIAD_COLOR, 'triad complete'), (~triad, NO_TRIAD_COLOR, 'triad unresolved')):
    centres, med, lo, hi = running_summary(mask)
    ax.fill_between(centres, lo, hi, color=color, alpha=0.20, linewidth=0, zorder=3)
    ax.plot(centres, med, color='white', linewidth=3.4, zorder=4, solid_capstyle='round')
    ax.plot(centres, med, color=color, linewidth=2.0, zorder=5, solid_capstyle='round', label=label)

ax.set_xlim(18, 100)   # no representative falls below 22.4% identity
ax.set_ylim(0, 103)
ax.set_xlabel('Sequence identity to best-matching reference (%)', fontsize=10.5, color=TEXT_DARK, labelpad=8)
ax.set_ylabel('Structural similarity to that reference (qTM-score, %)', fontsize=10.5, color=TEXT_DARK, labelpad=8)
# Set on two lines and centred in the band: the band is only 12 identity points
# wide, and one line of this length ran off the left edge of the axes.
ax.text((18 + TWILIGHT_MAX) / 2, 2.5, 'sequence\n"twilight zone"', ha='center', va='bottom', fontsize=8.6,
        color=TEXT_MUTED, style='italic', linespacing=1.45)
ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
for spine in ('left', 'bottom'):
    ax.spines[spine].set_color(RULE)
ax.tick_params(labelsize=9.5, colors=TEXT_MUTED)

# ---- right marginal: where the two groups actually separate
bin_edges = np.linspace(0, 103, 70)
for mask, color in ((~triad, NO_TRIAD_COLOR), (triad, TRIAD_COLOR)):
    density, _ = np.histogram(y[mask], bins=bin_edges, density=True)
    centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax_marg.fill_betweenx(centres, 0, density, color=color, alpha=0.30, linewidth=0)
    ax_marg.plot(density, centres, color=color, linewidth=1.4)
for threshold in (50, 85):
    ax_marg.axhline(threshold, color=RULE, linewidth=0.8, linestyle=(0, (4, 3)), zorder=1)
ax_marg.set_xlim(left=0)
ax_marg.set_xticks([])
ax_marg.tick_params(labelleft=False, length=0)
for spine in ax_marg.spines.values():
    spine.set_visible(False)
ax_marg.set_xlabel('density', fontsize=8.4, color=TEXT_MUTED, labelpad=6)

# ---- one legend, below the panel, outside the data
handles = [
    plt.Line2D([0], [0], color=TRIAD_COLOR, linewidth=2.0, marker='o', markersize=5, linestyle='-',
               label=f'Catalytic triad complete  ({n_triad:,})'),
    plt.Line2D([0], [0], color=NO_TRIAD_COLOR, linewidth=2.0, marker='o', markersize=5, linestyle='-',
               label=f'Triad not resolved by sequence alignment  ({n_no_triad:,})'),
    mpatches.Patch(facecolor='#B9C4CE', edgecolor='none', label='running median, IQR shaded'),
]
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.52, 0.004), ncol=3, frameon=False,
           fontsize=9.2, handlelength=2.0, columnspacing=2.6)

fig.savefig(OUT / 'seq_vs_structure.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'seq_vs_structure.pdf', facecolor='white')
print('saved', OUT / 'seq_vs_structure.png')
