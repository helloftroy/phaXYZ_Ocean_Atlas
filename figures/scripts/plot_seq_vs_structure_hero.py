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

A uniform-point-size variant was rendered once for comparison and dropped
per direct request -- bubble size (by genome count per paralog cluster) is
the version kept; render() still takes a use_size flag if that comparison
is ever wanted again.

Outputs:
    figures/seq_vs_structure_hero_sized.png / .pdf
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
# 3. sequence-identity-bin summary (the panel above the scatter): for
#    each bin, what fraction clears qTM>=50% and qTM>=85% -- answers
#    "as sequence identity falls, how often does structural support
#    remain?" directly, without having to read it off the scatter by eye
# ---------------------------------------------------------------------
BINS = [(0, 30, '<30%'), (30, 40, '30-40%'), (40, 50, '40-50%'), (50, 70, '50-70%'), (70, 100.0001, '≥70%')]
bin_stats = []
for lo, hi, label in BINS:
    sub = [r for r in rows if lo <= r['pident'] < hi]
    n = len(sub)
    p50 = 100 * sum(1 for r in sub if r['qtmscore'] >= 0.5) / n
    p85 = 100 * sum(1 for r in sub if r['qtmscore'] >= 0.85) / n
    bin_stats.append({'lo': lo, 'hi': min(hi, 100), 'label': label, 'n': n, 'p50': p50, 'p85': p85})
    print(f'  {label:8s} n={n:6,d}  qTM>=50%: {p50:5.1f}%   qTM>=85%: {p85:5.1f}%')

# ---------------------------------------------------------------------
# 4. the figure itself
# ---------------------------------------------------------------------

# Restyled 2026-09-29 per direct feedback that the original palette was not
# "drawing the eye in" -- checked, not just re-eyeballed: the original teal
# (#1E6E7A) FAILED the dataviz-skill palette validator's chroma-floor check
# outright ("reads gray", chroma 0.076 vs the >=0.1 floor), and the warm
# cream background (#FBFAF6) plus a warm-beige twilight-zone wash blended
# into an even muddier, lower-contrast combination on top of that. New teal
# (#0D9488) passes the validator cleanly against the orange it's paired
# with (chroma, CVD separation, normal-vision floor, contrast all PASS,
# node scripts/validate_palette.js "#0D9488,#C2622D" --mode light) and the
# background is now pure white, matching this project's other recent
# figures rather than the older cream-toned convention.
TRIAD_TRUE_COLOR = '#0D9488'
TRIAD_FALSE_COLOR = '#C2622D'
TRIAD_TRUE_LIGHT = '#7FCFC2'   # lighter tint of the same hue, for the ordinal bin-panel bars
BG = '#FFFFFF'
GRID_COLOR = '#E9ECEA'
TWILIGHT_ZONE_COLOR = '#EEF1F0'
BORDER_NEUTRAL = '#D5D9D7'
MUTED_TEXT = '#5B6E70'
DARK_TEXT = '#20302C'
n_triad_true = sum(1 for r in rows if r['triad_complete'])
n_triad_false = len(rows) - n_triad_true
max_n = max(r['n_genomes_in_cluster'] for r in rows)


def size_of(n, use_size):
    if not use_size:
        return 8
    # capped much smaller than a naive sqrt scale would give: the biggest
    # cluster (556 genomes) at a literal sqrt-proportional size drowned the
    # whole upper part of the plot in one solid blob (confirmed live on the
    # first render) -- area still grows with genome count, just compressed,
    # since legibility of the other 15,866 points matters more here than
    # exact linearity for the handful of largest clusters
    return 4 + 90 * (np.sqrt(n) / np.sqrt(max_n))


def render(use_size, out_stem):
    plt.rcParams.update({'font.family': 'DejaVu Sans'})
    fig = plt.figure(figsize=(15, 13), dpi=300)
    fig.patch.set_facecolor('white')
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 5.3], hspace=0.16, left=0.075, right=0.97, top=0.88, bottom=0.045)
    ax_bin = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1])

    # ---- bin panel -- evenly-spaced categorical groups, not spans matched
    # to actual bin width: the first version's variable-width blocks,
    # touching edge-to-edge with stark white divider lines, was confirmed
    # live to read as "big and ugly" -- real gaps between groups and a
    # narrower, more restrained bar width reads as a normal, polished
    # grouped bar chart instead
    ax_bin.set_facecolor(BG)
    n_bins = len(bin_stats)
    xpos = np.arange(n_bins)
    bar_w_outer, bar_w_inner = 0.52, 0.26
    ax_bin.grid(True, axis='y', color=GRID_COLOR, linewidth=1.2, zorder=0)
    ax_bin.set_axisbelow(True)
    for x, b in zip(xpos, bin_stats):
        ax_bin.bar(x, b['p50'], width=bar_w_outer, color=TRIAD_TRUE_LIGHT, zorder=2, edgecolor='none')
        ax_bin.bar(x, b['p85'], width=bar_w_inner, color=TRIAD_TRUE_COLOR, zorder=3, edgecolor='none')
        ax_bin.text(x, b['p50'] + 4, f"{b['p50']:.0f}%", ha='center', va='bottom', fontsize=9.5, color='#5B6E70')
        ax_bin.text(x, max(b['p85'] - 7, 4), f"{b['p85']:.0f}%", ha='center', va='bottom', fontsize=9.5,
                     color='white', fontweight='bold')
        ax_bin.text(x, -7, b['label'], ha='center', va='top', fontsize=10, color='#3A4A46', fontweight='bold')
        ax_bin.text(x, -17, f"n={b['n']:,}", ha='center', va='top', fontsize=8, color='#8B958F')
    ax_bin.set_xlim(-0.62, n_bins - 0.38)
    ax_bin.set_ylim(0, 112)
    ax_bin.set_xticks([])
    ax_bin.set_yticks([0, 50, 100])
    ax_bin.set_yticklabels(['0%', '50%', '100%'], fontsize=9)
    for spine in ('top', 'right', 'bottom'):
        ax_bin.spines[spine].set_visible(False)
    ax_bin.spines['left'].set_color(BORDER_NEUTRAL)
    ax_bin.tick_params(axis='y', length=0, colors='#8B958F')
    ax_bin.set_title('Structural support by sequence identity to best-matching reference', fontsize=12,
                      fontweight='bold', loc='left', color='#20302C', pad=10)
    legend_bin = [
        mpatches.Patch(color=TRIAD_TRUE_LIGHT, label='% with qTM ≥ 50% (structurally real)'),
        mpatches.Patch(color=TRIAD_TRUE_COLOR, label='% with qTM ≥ 85% (confident)'),
    ]
    ax_bin.legend(handles=legend_bin, loc='upper right', bbox_to_anchor=(1.0, 1.26), fontsize=9, frameon=False, ncol=2)

    # ---- main scatter ----
    ax.set_facecolor(BG)
    ax.axvspan(0, 30, color=TWILIGHT_ZONE_COLOR, alpha=1.0, zorder=0)
    ax.axhline(85, color='#B7BDB8', linewidth=1.0, linestyle=(0, (5, 3)), zorder=1)
    ax.axhline(50, color='#B7BDB8', linewidth=1.0, linestyle=(0, (5, 3)), zorder=1)

    # single combined scatter, ordered LARGEST-first (drawn on the bottom) so
    # small points aren't buried under big ones -- two separate per-color
    # calls (an earlier version of this figure) put every "triad complete"
    # point above every "incomplete" one regardless of size, which had the
    # same big-bubbles-hide-small-ones problem one layer up
    order = sorted(rows, key=lambda r: -r['n_genomes_in_cluster'])
    xs = [r['pident'] for r in order]
    ys = [100 * r['qtmscore'] for r in order]
    sizes = [size_of(r['n_genomes_in_cluster'], use_size) for r in order]
    colors = [TRIAD_TRUE_COLOR if r['triad_complete'] else TRIAD_FALSE_COLOR for r in order]
    ax.scatter(xs, ys, s=sizes, color=colors, alpha=0.4 if use_size else 0.35, linewidth=0, zorder=2)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 103)
    ax.set_xlabel('Sequence identity to best-matching reference (%)', fontsize=13, labelpad=10)
    ax.set_ylabel('Structural similarity to that reference (qTM-score, %)', fontsize=13, labelpad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, color=GRID_COLOR, linewidth=1.3, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=11)

    # reference-line labels back on the left -- the emptiest part of the
    # plot at these two heights (per direct request; an earlier version
    # moved them right to dodge a callout box that has since been removed)
    ref_label_kw = dict(fontsize=9, color='#5B6E70', ha='left', va='center', zorder=11,
                          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.82))
    ax.text(1.5, 85, 'qTM ≥ 85% (confident)', **ref_label_kw)
    ax.text(1.5, 50, 'qTM ≥ 50% (structurally real)', **ref_label_kw)
    ax.text(15, 101.5, 'sequence "twilight zone" (<30% identity)', fontsize=10.5, color='#7A8B86', ha='center',
            va='top', style='italic')

    # title / subtitle
    fig.text(0.075, 0.975, 'Sequence identity undersells it', fontsize=27, fontweight='bold', color='#20302C')
    fig.text(0.075, 0.952,
             f'{len(rows):,} representative phaC candidates, spanning nearly the whole verified dataset — {pct_confirmed:.1f}% fold as real phaC (qTM ≥ 50%)',
             fontsize=13.5, color='#5B6E70')
    fig.text(0.075, 0.936, 'even where sequence identity to any known reference is weak or absent.',
             fontsize=13.5, color='#5B6E70')

    # legends -- bottom-LEFT, under the twilight-zone band, the emptiest
    # part of the plot (confirmed live: bottom-right overlapped the point
    # cloud that spreads across the full x-range at low qTM)
    color_handles = [
        mpatches.Patch(color=TRIAD_TRUE_COLOR, label=f'Catalytic triad complete ({n_triad_true:,} / {len(rows):,})', alpha=0.75),
        mpatches.Patch(color=TRIAD_FALSE_COLOR, label=f'Triad not resolved by sequence alignment ({n_triad_false:,} / {len(rows):,})', alpha=0.75),
    ]
    leg_color = ax.legend(handles=color_handles, loc='lower left', bbox_to_anchor=(0.002, 0.002), fontsize=10,
                            frameon=True, facecolor='white', edgecolor=BORDER_NEUTRAL, title='Color', title_fontsize=10)
    ax.add_artist(leg_color)

    if use_size:
        size_vals = [1, 10, 100, max_n]
        size_handles = [mlines.Line2D([], [], marker='o', color='none', markerfacecolor='#5B6E70', markeredgecolor='#3A4A46',
                                       markeredgewidth=0.6, alpha=0.85, markersize=np.sqrt(size_of(v, True)), label=f'{v:,}')
                         for v in size_vals]
        ax.legend(handles=size_handles, loc='lower left', bbox_to_anchor=(0.002, 0.155), fontsize=9.5,
                  title='Genomes carrying this\nparalog cluster', title_fontsize=9.5, frameon=True,
                  facecolor='white', edgecolor=BORDER_NEUTRAL, labelspacing=1.3, borderpad=1.0, handletextpad=1.6)

    out_path = OUT / f'{out_stem}.png'
    fig.savefig(out_path, dpi=300, facecolor='white')
    fig.savefig(OUT / f'{out_stem}.pdf', facecolor='white')
    print('saved', out_path)
    plt.close(fig)


render(use_size=True, out_stem='seq_vs_structure_hero_sized')
