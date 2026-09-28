"""Hero summary figure for the comprehensive structural atlas (section 9.8):
one point per all_phac_dedup representative (15,867 of 15,880, QC-bad
targets dropped), x = sequence identity to its best-matching reference
(the same best_pident column used throughout this project), y =
structural similarity to that reference (qtmscore, the primary
structural-confidence metric established since section 9.7). Size = how
many genomes carry that representative's own phaC_cluster0.7 (70%-
identity) paralog cluster -- the same "how widespread is this lineage"
metric the cluster-landscape figure (section 5.6) uses, sqrt-scaled so
bubble AREA (not radius) is proportional. Color = catalytic-triad
completeness (catalytic_domain/phac_catalytic_triad.tsv's alignment-
column flag -- the only method with full 15,867-target coverage; the
more accurate structural geometric check from section 9.12 only applies
to the small subset with a folded PDB locally available, not this whole
comprehensive set, so is deliberately not mixed in here to keep one
consistent method across every point).

The plot earns the "hero" framing directly from its own data, not just
styling: qtmscore alone confirms 94.7% of the whole comprehensive set as
structurally real phaC (qtmscore>=0.5) even where sequence evidence was
weak -- and specifically, of the 1,043 representatives with <30%
sequence identity to their best reference (the "twilight zone" where
sequence-alone methods essentially can't tell homology from noise),
10.4% (108) are still structurally confident (qtmscore>=0.85) -- real
phaC that pure sequence search would have missed or discounted.

Usage:
    python figures/scripts/plot_seq_vs_structure_hero.py

Outputs:
    figures/seq_vs_structure_hero.png / .pdf
    figures/seq_vs_structure_hero.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.patheffects as pe
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

# ---------------------------------------------------------------------
# 1. load everything, one row per usable all_phac_dedup representative
# ---------------------------------------------------------------------
dedup_ids = set()
with open(ROOT / 'structure_prediction/all_phac_dedup_representatives.faa') as f:
    for line in f:
        if line.startswith('>'):
            dedup_ids.add(line[1:].split()[0])

bad_targets = _phac_qc.load_bad_targets()
targets = sorted(dedup_ids - bad_targets)
print(f'{len(dedup_ids):,} all_phac_dedup representatives -> {len(targets):,} after QC')

pident = {}
best_query = {}
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in dedup_ids:
            pident[row['target_id']] = float(row['best_pident'])
            best_query[row['target_id']] = row['best_query']

qtm = {}
with open(ROOT / 'structural_evidence_best_hit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['group'] == 'all_phac_dedup':
            qtm[row['candidate_id']] = float(row['qtmscore'])

triad = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in dedup_ids:
            triad[row['target_id']] = row['triad_complete'] == 'True'

# cluster0.7 assignment + genome count per cluster (same method as
# figures/scripts/plot_phac_cluster_landscape.py, section 5.6)
assignments = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if row[1] in dedup_ids:
            assignments[row[1]] = row[0]

cluster_of_interest = set(assignments.values())
print('building full cluster0.7 -> genome-count map (reused across all representatives) ...')
target_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        target_to_cluster[row[1]] = row[0]

genome_of_target = {}
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        genome_of_target[row['target_id']] = row['genome']

cluster_genome_set = defaultdict(set)
for t, c in target_to_cluster.items():
    if c in cluster_of_interest and t not in bad_targets:
        g = genome_of_target.get(t)
        if g:
            cluster_genome_set[c].add(g)
n_genomes_of_cluster = {c: len(gs) for c, gs in cluster_genome_set.items()}

# ---------------------------------------------------------------------
# 2. assemble the plotting table
# ---------------------------------------------------------------------
rows = []
for t in targets:
    if t not in pident or t not in qtm:
        continue
    c = assignments.get(t)
    n_g = n_genomes_of_cluster.get(c, 1)
    rows.append({
        'target_id': t, 'pident': pident[t], 'qtmscore': qtm[t],
        'triad_complete': triad.get(t, False), 'cluster_id': c, 'n_genomes_in_cluster': n_g,
    })
print(f'{len(rows):,} points with full data')

out_tsv = OUT / 'seq_vs_structure_hero.tsv'
with open(out_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print('wrote', out_tsv)

# twilight-zone headline stat
low_seq = [r for r in rows if r['pident'] < 30]
low_seq_confident = [r for r in low_seq if r['qtmscore'] >= 0.85]
print(f'twilight zone (<30% seq identity): {len(low_seq)}, of which structurally confident (qtm>=0.85): {len(low_seq_confident)} '
      f'({100*len(low_seq_confident)/len(low_seq):.1f}%)')
pct_confirmed = 100 * sum(1 for r in rows if r['qtmscore'] >= 0.5) / len(rows)
print(f'overall qtmscore>=0.5: {pct_confirmed:.1f}%')

# ---------------------------------------------------------------------
# 3. the figure itself
# ---------------------------------------------------------------------
TRIAD_TRUE_COLOR = '#1E6E7A'
TRIAD_FALSE_COLOR = '#C2622D'
BG = '#FBFAF6'

plt.rcParams.update({'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(15, 11.5), dpi=300)
fig.patch.set_facecolor('white')
ax.set_facecolor(BG)

# twilight-zone shaded band, drawn first (behind everything)
ax.axvspan(0, 30, color='#EDE6D6', alpha=0.55, zorder=0)
ax.axhline(85, color='#B7BDB8', linewidth=1.0, linestyle=(0, (5, 3)), zorder=1)
ax.axhline(50, color='#B7BDB8', linewidth=1.0, linestyle=(0, (5, 3)), zorder=1)

max_n = max(r['n_genomes_in_cluster'] for r in rows)


def size_of(n):
    # capped much smaller than a naive sqrt scale would give: the biggest
    # cluster (556 genomes) at a literal sqrt-proportional size drowned the
    # whole upper part of the plot in one solid blob (confirmed live on the
    # first render) -- area still grows with genome count, just compressed,
    # since legibility of the other 15,866 points matters more here than
    # exact linearity for the handful of largest clusters
    return 4 + 90 * (np.sqrt(n) / np.sqrt(max_n))


# single combined scatter, ordered LARGEST-first (drawn on the bottom) so
# small points aren't buried under big ones -- two separate per-color calls
# (the first version of this figure) put every "triad complete" point above
# every "incomplete" one regardless of size, which had the same
# big-bubbles-hide-small-ones problem one layer up
order = sorted(rows, key=lambda r: -r['n_genomes_in_cluster'])
xs = [r['pident'] for r in order]
ys = [100 * r['qtmscore'] for r in order]
sizes = [size_of(r['n_genomes_in_cluster']) for r in order]
colors = [TRIAD_TRUE_COLOR if r['triad_complete'] else TRIAD_FALSE_COLOR for r in order]
ax.scatter(xs, ys, s=sizes, color=colors, alpha=0.4, linewidth=0, zorder=2)

ax.set_xlim(0, 100)
ax.set_ylim(0, 103)
ax.set_xlabel('Sequence identity to best-matching reference (%)', fontsize=13, labelpad=10)
ax.set_ylabel('Structural similarity to that reference (qTM-score, %)', fontsize=13, labelpad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, color='white', linewidth=1.3, zorder=0)
ax.set_axisbelow(True)
ax.tick_params(labelsize=11)

# reference-line labels sit on the LEFT (the twilight-zone band is the
# emptiest part of the plot at high qTM) so they don't compete with the
# legends, which live bottom-right where the point cloud thins out
ax.text(1, 86.3, 'qTM ≥ 85% (confident)', fontsize=9, color='#7A8580', ha='left', va='bottom')
ax.text(1, 51.3, 'qTM ≥ 50% (structurally real)', fontsize=9, color='#7A8580', ha='left', va='bottom')
ax.text(15, 101.5, 'sequence "twilight zone" (<30% identity)', fontsize=10.5, color='#9A8A60', ha='center', va='top',
        style='italic')

# twilight-zone callout: placed in the one clearly empty patch of the plot
# (low sequence identity, low-to-mid structural similarity -- few real
# candidates land there since low qTM there would mean neither sequence
# nor structure supports the call) rather than on top of the dense cloud
ax.annotate(f'{len(low_seq_confident)} of {len(low_seq):,} twilight-zone candidates ({100*len(low_seq_confident)/len(low_seq):.0f}%)\nare still structurally confident (qTM ≥ 85%) --\nreal phaC that sequence search alone would miss.',
            xy=(27, 88), xytext=(20, 14), fontsize=11.5, color='#3A3226', ha='left', va='bottom',
            arrowprops=dict(arrowstyle='-|>', color='#5B6E70', lw=1.4, connectionstyle='arc3,rad=0.25'),
            bbox=dict(boxstyle='round,pad=0.55', facecolor='white', edgecolor='#D8D2C0', alpha=0.96), zorder=10)

# title / subtitle
fig.text(0.085, 0.965, 'Sequence identity undersells it', fontsize=27, fontweight='bold', color='#20302C')
fig.text(0.085, 0.935,
         f'{len(rows):,} representative phaC candidates, spanning nearly the whole verified dataset — {pct_confirmed:.1f}% fold as real phaC (qTM ≥ 50%)',
         fontsize=13.5, color='#5B6E70')
fig.text(0.085, 0.915, 'even where sequence identity to any known reference is weak or absent.',
         fontsize=13.5, color='#5B6E70')

# legends
size_vals = [1, 10, 100, max_n]
size_handles = [mlines.Line2D([], [], marker='o', color='none', markerfacecolor='#5B6E70', markeredgecolor='#3A4A46',
                               markeredgewidth=0.6, alpha=0.85, markersize=np.sqrt(size_of(v)), label=f'{v:,}')
                 for v in size_vals]
leg_size = ax.legend(handles=size_handles, loc='lower right', bbox_to_anchor=(0.998, 0.14), fontsize=9.5,
                      title='Genomes carrying this\nparalog cluster', title_fontsize=9.5, frameon=True,
                      facecolor='white', edgecolor='#D8D2C0', labelspacing=1.3, borderpad=1.0, handletextpad=1.6)
ax.add_artist(leg_size)

color_handles = [
    mpatches.Patch(color=TRIAD_TRUE_COLOR, label='Catalytic triad complete', alpha=0.75),
    mpatches.Patch(color=TRIAD_FALSE_COLOR, label='Triad incomplete (alignment-column check)', alpha=0.75),
]
ax.legend(handles=color_handles, loc='lower right', bbox_to_anchor=(0.998, 0.005), fontsize=10, frameon=True,
          facecolor='white', edgecolor='#D8D2C0', title='Color', title_fontsize=10)

fig.text(0.085, 0.018,
         'One point per all_phac_dedup representative (near-duplicate-collapsed; large diverse paralog clusters contribute more than one representative). Triad status is the\n'
         "alignment-column check (full coverage) -- section 9.12's own structural geometric check is more accurate but only applies where a folded PDB is locally available.",
         fontsize=8.3, color='#8B958F')

fig.subplots_adjust(left=0.075, right=0.97, top=0.885, bottom=0.095)

out_path = OUT / 'seq_vs_structure_hero.png'
fig.savefig(out_path, dpi=300, facecolor='white')
fig.savefig(OUT / 'seq_vs_structure_hero.pdf', facecolor='white')
print('saved', out_path)
