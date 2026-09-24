"""Summary scatter for catalytic_domain/audit_divergent35_structural_triad.py's
output: for every figures/phac_divergent35_validated.tsv candidate with a
folded PDB, plots the two key catalytic-triad hydrogen-bond distances
(nucleophile-to-His, His-to-Asp) against each other, colored by
confidence tier and shaped by which nucleophile the best geometric match
used. See figures/PHA_CLEAN_RESULTS.md for the full write-up -- the
headline: independent HMM/pathway confidence tier tracks structural
triad plausibility closely (92% of HIGH-tier candidates land in the
tight corner vs. 20-25% of MEDIUM/LOW), and a real minority of the
no-confirmed-Cys candidates use a genuine Ser (rarely Thr) nucleophile
instead, not a Cys the HMM approach simply missed.

Usage:
    python figures/scripts/plot_divergent35_triad_geometry.py

Output:
    figures/divergent35_triad_geometry.png/.pdf
"""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
TIGHT_CUTOFF = 4.5

TIER_COLOR = {'HIGH': '#1E6E7A', 'MEDIUM': '#C9A227', 'LOW': '#9E3B3B'}
NUC_MARKER = {'Cys': 'o', 'S': '^', 'T': 's', 'none': 'x'}
NUC_LABEL = {'Cys': 'Cys (canonical)', 'S': 'Ser (alternative)', 'T': 'Thr (alternative)', 'none': 'no Ser/Thr/Cys at all'}

rows = []
with open(OUT / 'divergent35_structural_triad_audit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        rows.append(row)

fig, ax = plt.subplots(figsize=(9, 7.5), dpi=300)

max_val = 35
tight_box = mpatches.Rectangle((0, 0), TIGHT_CUTOFF, TIGHT_CUTOFF, facecolor='#DCEDEA', edgecolor='none', zorder=0)
ax.add_patch(tight_box)
ax.annotate(f'"tight" region\n(<= {TIGHT_CUTOFF}A both legs)', xy=(TIGHT_CUTOFF, TIGHT_CUTOFF), xytext=(15, 12),
            fontsize=9, color='#3A6B63', ha='left',
            arrowprops=dict(arrowstyle='->', color='#3A6B63', lw=1.2))

for row in rows:
    if row['nucleophile'] == 'none':
        continue
    x = float(row['nuc_his_dist_A'])
    y = float(row['his_asp_dist_A'])
    tier = row['confidence_tier'].split(':')[0]
    marker = NUC_MARKER.get(row['nucleophile'], 'x')
    ax.scatter(min(x, max_val), min(y, max_val), s=70, marker=marker, color=TIER_COLOR.get(tier, '#888'),
               edgecolor='#20302C', linewidth=0.6, alpha=0.85, zorder=3)

ax.set_xlim(0, max_val)
ax.set_ylim(0, max_val)
ax.set_xlabel('nucleophile <-> His distance (A)  [closest SG/OG/OG1 to closest NE2/ND1]')
ax.set_ylabel('His <-> Asp distance (A)  [closest NE2/ND1 to closest OD1/OD2]')
ax.set_title('Structural catalytic-triad geometry across the divergent35 set\n(most sequence-divergent validated phaC candidates, <35% identity to nearest ref.)',
              fontsize=10.5, fontweight='bold')

tier_handles = [mpatches.Patch(color=c, label=f'{t} confidence') for t, c in TIER_COLOR.items()]
nuc_handles = [mlines.Line2D([], [], color='#555', marker=m, linestyle='None', markersize=8, label=NUC_LABEL[n])
               for n, m in NUC_MARKER.items() if n != 'none']
leg1 = ax.legend(handles=tier_handles, loc='upper right', title='HMM/pathway confidence tier', fontsize=8, title_fontsize=8.5, frameon=False)
ax.add_artist(leg1)
ax.legend(handles=nuc_handles, loc='lower right', title='best-fit nucleophile', fontsize=8, title_fontsize=8.5, frameon=False)

fig.tight_layout()
fig.savefig(OUT / 'divergent35_triad_geometry.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'divergent35_triad_geometry.pdf', facecolor='white')
print('saved figures/divergent35_triad_geometry.png/.pdf')
