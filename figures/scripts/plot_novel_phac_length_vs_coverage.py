"""Are the 255 novel phaC candidates compact proteins, or truncated gene calls?

Protein length against the fraction of its best structural reference that the
Foldseek alignment covers (ESMFold models; PHA_CLEAN_RESULTS.md section 13.4).
The two readings separate cleanly on this plot:

  truncated gene call   a fragment aligns well over only part of its reference,
                        so it sits low on the y axis however long it is.
  compact architecture  a short protein that still spans its whole reference
                        sits at the top left, which is where this set is.

Reference coverage rather than query coverage on the y axis, deliberately: the
query side cannot distinguish the two cases. A fragment aligns over most of
*itself* by construction, so qcov is high for a truncation; only the reference
side shows that half the fold is missing.

Colour is the TM-score over the alignment, which is the fold-match strength, so
the panel carries all three quantities the section's claim rests on.

Usage:
    python figures/scripts/plot_novel_phac_length_vs_coverage.py
"""
import csv
import statistics
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
EVIDENCE = ROOT / 'structural_evidence_best_hit.tsv'
ANNOTATIONS = ROOT / 'phac_recovery/novel_phac_new_candidates.tsv'

COMPLETE_CUTOFF = 0.75
TEXT_DARK, TEXT_MUTED, RULE, GRID = '#1E2630', '#5B6670', '#B8C0C8', '#EAEDF1'

annotations = {r['gene_id']: r for r in csv.DictReader(open(ANNOTATIONS, newline=''), delimiter='\t')}
rows = []
with open(EVIDENCE, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['group'] != 'novel_phac_candidates' or row['status'] != 'hit':
            continue
        annotation = annotations.get(row['candidate_id'])
        if annotation is None:
            continue
        rows.append(dict(gene_id=row['candidate_id'], length=int(annotation['gene_len']),
                         tcov=float(row['tcov_struct']), tm=float(row['alntmscore']),
                         pident=float(annotation['pident'])))
print(f'{len(rows)} candidates')

length = np.array([r['length'] for r in rows])
tcov = np.array([r['tcov'] for r in rows])
tm = np.array([r['tm'] for r in rows])
median_length = statistics.median(length)
median_tcov = statistics.median(tcov)
n_complete = int((tcov >= COMPLETE_CUTOFF).sum())
print(f'median length {median_length:.0f} aa, median reference coverage {100 * median_tcov:.1f}%, '
      f'{n_complete}/{len(rows)} at or above {100 * COMPLETE_CUTOFF:.0f}%')

plt.rcParams.update({'font.family': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 10})
fig, ax = plt.subplots(figsize=(8.0, 6.0), dpi=300)

ax.axhline(COMPLETE_CUTOFF, color=RULE, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
ax.axhline(median_tcov, color='#0d366b', linewidth=1.0, linestyle=(0, (1, 2)), zorder=1)
ax.axvline(median_length, color='#0d366b', linewidth=1.0, linestyle=(0, (1, 2)), zorder=1)

scatter = ax.scatter(length, tcov, c=tm, cmap='YlGnBu', vmin=0.80, vmax=0.98,
                     s=42, edgecolors='white', linewidths=0.6, zorder=3)

ax.set_xlabel('Candidate synthase length (aa)', fontsize=11, color=TEXT_DARK, labelpad=8)
ax.set_ylabel('Fraction of best reference covered by the structural alignment',
              fontsize=11, color=TEXT_DARK, labelpad=8)
ax.set_xlim(220, 620)
ax.set_ylim(0.25, 1.04)
ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color(RULE)
ax.tick_params(colors=TEXT_MUTED, length=3, labelsize=9.5)

ax.text(624, COMPLETE_CUTOFF, f'{COMPLETE_CUTOFF:.2f}', ha='left', va='center', fontsize=8.6, color=TEXT_MUTED)
ax.text(median_length, 1.045, f'median {median_length:.0f} aa', ha='center', va='bottom',
        fontsize=9, color='#0d366b')
ax.text(228, median_tcov + 0.012, f'median {100 * median_tcov:.1f}%', ha='left', va='bottom',
        fontsize=9, color='#0d366b')
ax.text(228, 0.29,
        f'{n_complete} of {len(rows)} cover at least {100 * COMPLETE_CUTOFF:.0f}% of their reference.\n'
        'A truncated gene call would sit along the bottom.',
        ha='left', va='bottom', fontsize=9.2, color=TEXT_DARK, linespacing=1.5)

cbar = fig.colorbar(scatter, ax=ax, pad=0.09, fraction=0.045)
cbar.set_label('TM-score over the structural alignment', fontsize=9.5, color=TEXT_DARK)
cbar.ax.tick_params(labelsize=8.5, colors=TEXT_MUTED, length=2)
cbar.outline.set_visible(False)

fig.tight_layout()
fig.savefig(OUT / 'novel_phac_length_vs_coverage.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'novel_phac_length_vs_coverage.pdf', facecolor='white')
print('saved', OUT / 'novel_phac_length_vs_coverage.png')
