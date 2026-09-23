"""%identity vs. query coverage to the audited phaC reference set (1,723
confirmed-genuine references as of the 2026-09-22 fix -- was 1,875, but
all 152 newly-excluded accessions from that fix were themselves part of
this reference set; see figures/PHA_CLEAN_RESULTS.md section 2), one
point per verified phaC target, colored by the five evidence tiers from
section 3 (catalytic triad / HMM / pathway-richness).

Motivation: a high %identity alone is not strong evidence if it's only
over a tiny fraction of the protein -- e.g. 50% identity over 8% of the
candidate's length is a short, possibly coincidental, conserved-motif-only
match, not evidence of a full-length homolog. Plotting pident against
qcov (not pident alone) makes that distinction visible directly.
"""
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'

ORDER = ['Catalytic triad complete', 'HMM-supported (no triad)', 'No HMM/triad, >=5 other PHA genes',
         'No HMM/triad, 1-4 other PHA genes', 'phaC only']
LABELS = {
    'Catalytic triad complete': 'Catalytic triad complete',
    'HMM-supported (no triad)': 'HMM-supported (no triad)',
    'No HMM/triad, >=5 other PHA genes': '≥5 other PHA genes',
    'No HMM/triad, 1-4 other PHA genes': '1-4 other PHA genes',
    'phaC only': 'phaC only',
}
COLORS = {
    'Catalytic triad complete': '#1E6E7A',
    'HMM-supported (no triad)': '#3E8914',
    'No HMM/triad, >=5 other PHA genes': '#7A9B3E',
    'No HMM/triad, 1-4 other PHA genes': '#C2A83E',
    'phaC only': '#9E3B3B',
}
ZORDER = {k: 5 - i for i, k in enumerate(ORDER)}  # plot strongest tier on top

points = defaultdict(list)  # tier -> list of (qcov_pct, pident)
n_no_hit = 0
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        if not row['pident']:
            n_no_hit += 1
            continue
        tier = row['evidence_tier']
        points[tier].append((100 * float(row['qcov']), float(row['pident'])))

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, ax = plt.subplots(figsize=(11.5, 8), dpi=300)

for tier in sorted(ORDER, key=lambda t: ZORDER[t], reverse=True):
    pts = points[tier]
    if not pts:
        continue
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.scatter(xs, ys, s=7, alpha=0.35, color=COLORS[tier], edgecolors='none',
               zorder=ZORDER[tier], label=f'{LABELS[tier]} (n={len(pts):,})')

# reference marker for the "weak evidence" example named in the request
ax.scatter([8], [50], marker='x', s=140, color='black', zorder=10, linewidths=2.5)
ax.annotate('50% identity,\n8% coverage', xy=(8, 50), xytext=(18, 38),
            fontsize=9, ha='left', arrowprops=dict(arrowstyle='->', color='black', lw=1))

ax.set_xlabel('Query coverage to best-matching audited reference (%)', fontsize=12)
ax.set_ylabel('% identity to best-matching audited reference', fontsize=12)
ax.set_xlim(-2, 102)
ax.set_ylim(0, 102)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(color='#E4E7E2', linewidth=0.6, zorder=0)
ax.set_axisbelow(True)
leg = ax.legend(loc='upper left', fontsize=9.5, frameon=False, markerscale=3)

fig.suptitle('phaC Sequence Evidence: Identity vs. Coverage to the Audited Reference Set', fontsize=14.5, fontweight='bold', y=0.985)
fig.text(0.5, 0.945,
          f'{sum(len(v) for v in points.values()):,} of {sum(len(v) for v in points.values()) + n_no_hit:,} verified phaC targets have >=1 hit (e<1e-5) to the 1,723-sequence audited reference set; the rest ({n_no_hit:,}) have none at all.',
          ha='center', fontsize=9.3, color='#5B6E70')

footnote = (
    'Reference set = 1,723 phaC/FUSION references confirmed genuine by the full InterPro domain audit (PHA_CLEAN_RESULTS.md §2), not just the\n'
    '14-sequence "trusted" benchmark used earlier. High identity concentrated at low coverage (bottom-right-ish) marks short, possibly-coincidental\n'
    'conserved-motif-only matches rather than full-length homology -- coverage matters as much as identity for judging real evidence strength.'
)
fig.text(0.5, 0.01, footnote, ha='center', va='bottom', fontsize=8.0, color='#5B6E70')

fig.subplots_adjust(left=0.09, right=0.97, top=0.89, bottom=0.13)

out_path = OUT / 'phac_pident_vs_coverage.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_pident_vs_coverage.pdf', facecolor='white')
