"""How much phaC cluster diversity is left to find, and how much the answer
depends on which studies are in the pool.

Replaces the phylum-point overlay (plot_phac_diversity_phylum_landscape.py).
The per-phylum points were hard to read at a glance and answered a different
question -- which lineages are rich -- than the curve itself does. What is
left is the saturation question, asked with two sampling units instead of one:

  genome order   genomes drawn one at a time in random order. The standard
                 rarefaction curve, and the optimistic one: consecutive draws
                 routinely come from the same study, so the curve benefits
                 from within-study redundancy that a genuinely new sample
                 would not provide.
  study order    whole studies drawn in random order, every genome in a study
                 entering at once, plotted against the same cumulative-genome
                 axis. Genomes within a study share a cruise, a depth range, a
                 protocol and an assembly pipeline, so they are not
                 independent draws from the ocean. This treats the study as
                 the unit of independent sampling, which is the conservative
                 reading.

The gap between them is the figure's point: it is how much of the apparent
discovery rate is an artifact of pooling studies rather than a property of
phaC diversity. Where the study curve sits below the genome curve at the same
number of genomes, those genomes were buying less new diversity than a random
draw implies.

Chao2 (Chao 1987, incidence-based, genome as the incidence unit) is drawn as a
horizontal asymptote with a log-normal 95% interval, against the observed
count. Chao2 is also reported at the other mmseqs identity thresholds this
project clustered at (50%, 60%, 90%), since the saturation fraction depends on
where the cluster boundary is drawn and quoting it at 70% alone would overstate
how threshold-independent the result is.

Usage:
    python figures/scripts/plot_phac_accumulation_curves.py
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

RANDOM_SEED = 42
N_PERM = 300
N_CURVE_POINTS = 240
THRESHOLDS = ['0.5', '0.6', '0.7', '0.9']
PRIMARY = '0.7'

GENOME_COLOR, GENOME_BAND = '#0D9488', '#CCE8E4'
STUDY_COLOR, STUDY_BAND = '#C2622D', '#F3DCCB'
CHAO_COLOR, CHAO_BAND = '#3A4A46', '#E8EAE6'
TEXT_DARK, TEXT_MUTED, GRID = '#20302C', '#5B6E70', '#EBEEEC'

rng = np.random.default_rng(RANDOM_SEED)
bad_targets = _phac_qc.load_bad_targets()

target_to_genomes = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        if target_id not in bad_targets:
            target_to_genomes[target_id].add(genome)


def genome_clusters_at(threshold):
    """genome -> set of cluster ids, at one mmseqs identity threshold."""
    mapping = defaultdict(set)
    with open(FA / f'phaC_cluster{threshold}_cluster.tsv', newline='') as f:
        for cluster_id, target_id in csv.reader(f, delimiter='\t'):
            if target_id in bad_targets:
                continue
            for genome in target_to_genomes.get(target_id, ()):
                mapping[genome].add(cluster_id)
    return mapping


def chao2(cluster_incidence, n_units):
    """Chao2 with a log-normal 95% interval (Chao 1987; Colwell's EstimateS
    formulation). cluster_incidence maps cluster -> number of sampling units it
    was found in. Q1/Q2 are the clusters found in exactly one and exactly two."""
    observed = len(cluster_incidence)
    q1 = sum(1 for c in cluster_incidence.values() if c == 1)
    q2 = sum(1 for c in cluster_incidence.values() if c == 2)
    a = (n_units - 1) / n_units
    if q2 > 0:
        estimate = observed + a * q1 ** 2 / (2 * q2)
        ratio = q1 / q2
        variance = q2 * (a / 2 * ratio ** 2 + a ** 2 * ratio ** 3 + a ** 2 / 4 * ratio ** 4)
    else:
        # Bias-corrected form, used when no cluster is seen in exactly two units.
        estimate = observed + a * q1 * (q1 - 1) / 2
        variance = a * q1 * (q1 - 1) / 2 + a * q1 * (2 * q1 - 1) ** 2 / 4 - a * q1 ** 4 / (4 * estimate)
    extra = estimate - observed
    if extra <= 0 or variance <= 0:
        return estimate, observed, estimate, q1, q2
    k = math.exp(1.96 * math.sqrt(math.log(1 + variance / extra ** 2)))
    return estimate, observed + extra / k, observed + extra * k, q1, q2


# ------------------------------------------------------------------ primary set
genome_clusters = genome_clusters_at(PRIMARY)
genomes = sorted(genome_clusters)
cluster_sets = [genome_clusters[g] for g in genomes]
n_genomes = len(genomes)
n_clusters = len(set().union(*cluster_sets))

study_of = {}
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        study_of[row['genome']] = row['study_id']
study_members = defaultdict(list)
for i, genome in enumerate(genomes):
    study_members[study_of.get(genome, 'unknown')].append(i)
studies = sorted(study_members)
print(f'{n_genomes:,} phaC-positive genomes, {n_clusters:,} clusters at {PRIMARY}, {len(studies)} studies')
print(f'largest study {max(len(v) for v in study_members.values()):,} genomes, '
      f'median {int(np.median([len(v) for v in study_members.values()])):,}')

xs = np.unique(np.round(np.linspace(1, n_genomes, N_CURVE_POINTS)).astype(int))

# genome-order rarefaction
genome_curves = np.zeros((N_PERM, len(xs)))
order = np.arange(n_genomes)
for p in range(N_PERM):
    seen, xi = set(), 0
    for i, gi in enumerate(rng.permutation(order), start=1):
        seen.update(cluster_sets[gi])
        while xi < len(xs) and i == xs[xi]:
            genome_curves[p, xi] = len(seen)
            xi += 1

# study-block accumulation, on the same cumulative-genome axis. Each permutation
# gives a step function (a whole study lands at once); it is read onto the common
# grid by step lookup, not linear interpolation, so the curve is not smoothed
# across a block boundary that the data says is a jump.
study_curves = np.zeros((N_PERM, len(xs)))
study_index = np.arange(len(studies))
for p in range(N_PERM):
    seen = set()
    cumulative_genomes, cumulative_clusters = [0], [0]
    total = 0
    for si in rng.permutation(study_index):
        members = study_members[studies[si]]
        for i in members:
            seen.update(cluster_sets[i])
        total += len(members)
        cumulative_genomes.append(total)
        cumulative_clusters.append(len(seen))
    positions = np.searchsorted(np.array(cumulative_genomes), xs, side='right') - 1
    study_curves[p] = np.array(cumulative_clusters)[positions]

genome_mean = genome_curves.mean(axis=0)
genome_lo, genome_hi = np.percentile(genome_curves, [2.5, 97.5], axis=0)
study_mean = study_curves.mean(axis=0)
study_lo, study_hi = np.percentile(study_curves, [2.5, 97.5], axis=0)

half = len(xs) // 2
print(f'at {xs[half]:,} genomes: genome order {genome_mean[half]:,.0f} clusters, '
      f'study order {study_mean[half]:,.0f} ({100 * study_mean[half] / genome_mean[half]:.0f}% of it)')

# ------------------------------------------------------------------ Chao2
incidence = defaultdict(int)
for s in cluster_sets:
    for c in s:
        incidence[c] += 1
estimate, ci_lo, ci_hi, q1, q2 = chao2(incidence, n_genomes)
print(f'Chao2 at {PRIMARY}: {estimate:,.0f} (95% CI {ci_lo:,.0f}-{ci_hi:,.0f}); '
      f'observed {n_clusters:,} = {100 * n_clusters / estimate:.0f}%; Q1={q1:,} Q2={q2:,}')

rows = []
for threshold in THRESHOLDS:
    mapping = genome_clusters_at(threshold) if threshold != PRIMARY else genome_clusters
    counts = defaultdict(int)
    for s in mapping.values():
        for c in s:
            counts[c] += 1
    est, lo, hi, t_q1, t_q2 = chao2(counts, len(mapping))
    rows.append(dict(threshold=threshold, n_genomes=len(mapping), observed=len(counts), chao2=est,
                     ci_lo=lo, ci_hi=hi, pct=100 * len(counts) / est, q1=t_q1, q2=t_q2))
    print(f'  {int(float(threshold) * 100)}% identity: observed {len(counts):,}, Chao2 {est:,.0f} '
          f'({lo:,.0f}-{hi:,.0f}), {100 * len(counts) / est:.0f}% sampled')

with open(OUT / 'phac_accumulation_chao2_by_threshold.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['identity_threshold', 'n_genomes', 'observed_clusters', 'chao2', 'chao2_ci_lo', 'chao2_ci_hi',
                'pct_observed', 'Q1_singletons', 'Q2_doubletons'])
    for row in rows:
        w.writerow([row['threshold'], row['n_genomes'], row['observed'], f"{row['chao2']:.0f}",
                    f"{row['ci_lo']:.0f}", f"{row['ci_hi']:.0f}", f"{row['pct']:.1f}", row['q1'], row['q2']])

with open(OUT / 'phac_accumulation_curves.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['genomes_sampled', 'genome_order_mean', 'genome_order_lo', 'genome_order_hi',
                'study_order_mean', 'study_order_lo', 'study_order_hi'])
    for i, x in enumerate(xs):
        w.writerow([x, f'{genome_mean[i]:.1f}', f'{genome_lo[i]:.1f}', f'{genome_hi[i]:.1f}',
                    f'{study_mean[i]:.1f}', f'{study_lo[i]:.1f}', f'{study_hi[i]:.1f}'])

# ------------------------------------------------------------------ figure
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5})
fig, ax = plt.subplots(figsize=(12.6, 8.2), dpi=300)

ax.axhspan(ci_lo, ci_hi, color=CHAO_BAND, zorder=0, linewidth=0)
ax.axhline(estimate, color=CHAO_COLOR, linewidth=1.5, linestyle=(0, (5, 3)), zorder=3)

ax.fill_between(xs, study_lo, study_hi, color=STUDY_BAND, zorder=1, linewidth=0)
ax.fill_between(xs, genome_lo, genome_hi, color=GENOME_BAND, zorder=1, linewidth=0)
ax.plot(xs, study_mean, color=STUDY_COLOR, linewidth=2.4, zorder=4)
ax.plot(xs, genome_mean, color=GENOME_COLOR, linewidth=2.4, zorder=4)

ax.scatter([n_genomes], [n_clusters], s=90, color=GENOME_COLOR, edgecolors='white', linewidths=1.4, zorder=6)

ax.set_xlim(0, n_genomes * 1.19)
ax.set_ylim(0, ci_hi * 1.06)
ax.set_xlabel('phaC-positive genomes sampled', fontsize=12, color=TEXT_DARK)
ax.set_ylabel(f'Distinct phaC clusters recovered (70% identity)', fontsize=12, color=TEXT_DARK)
ax.grid(True, color=GRID, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
for side in ('top', 'right'):
    ax.spines[side].set_visible(False)
for side in ('left', 'bottom'):
    ax.spines[side].set_color('#C3C2B7')
ax.xaxis.set_major_formatter(lambda v, _: f'{v:,.0f}')

ax.annotate(f'Chao2 estimate {estimate:,.0f}\n95% CI {ci_lo:,.0f}–{ci_hi:,.0f}',
            xy=(n_genomes * 1.01, estimate), xytext=(n_genomes * 1.02, estimate),
            ha='left', va='center', fontsize=9.8, color=CHAO_COLOR)
ax.annotate(f'{n_clusters:,} observed\n{100 * n_clusters / estimate:.0f}% of Chao2',
            xy=(n_genomes, n_clusters), xytext=(n_genomes * 1.02, n_clusters),
            ha='left', va='center', fontsize=9.8, color=GENOME_COLOR, fontweight='bold')

# Mark the gap mid-curve, not at the end: both curves necessarily meet at the full
# sample (the last study drawn contributes every cluster still missing), so the
# endpoint carries no information and labelling it would invite the reading that
# the two sampling units agree.
gap_x = xs[half]
ax.plot([gap_x, gap_x], [study_mean[half], genome_mean[half]], color=TEXT_MUTED, linewidth=1.0,
        linestyle=(0, (2, 2)), zorder=5)
for y in (study_mean[half], genome_mean[half]):
    ax.plot([gap_x - n_genomes * 0.004, gap_x + n_genomes * 0.004], [y, y], color=TEXT_MUTED,
            linewidth=1.0, zorder=5)
ax.annotate(f'At {gap_x:,} genomes the study-order curve\nhas found {100 * study_mean[half] / genome_mean[half]:.0f}% '
            f'of what genome order claims\n({study_mean[half]:,.0f} vs {genome_mean[half]:,.0f} clusters)',
            xy=(gap_x, (study_mean[half] + genome_mean[half]) / 2),
            xytext=(gap_x - n_genomes * 0.03, (study_mean[half] + genome_mean[half]) / 2 + 620),
            ha='right', va='center', fontsize=9.3, color=TEXT_DARK, linespacing=1.45,
            arrowprops=dict(arrowstyle='-', color=TEXT_MUTED, linewidth=0.9,
                            connectionstyle='angle3,angleA=0,angleB=70'))

handles = [
    plt.Line2D([0], [0], color=GENOME_COLOR, linewidth=2.4,
               label='Sampling unit: genome\n(genomes in random order)'),
    plt.Line2D([0], [0], color=STUDY_COLOR, linewidth=2.4,
               label=f'Sampling unit: study\n(all {len(studies)} studies in random order,\nevery genome entering at once)'),
    Patch(facecolor='#D5DED9', edgecolor='none', label=f'95% interval over {N_PERM} permutations'),
    plt.Line2D([0], [0], color=CHAO_COLOR, linewidth=1.5, linestyle=(0, (5, 3)),
               label='Chao2 asymptotic richness (95% CI shaded)'),
]
# Anchored below the Chao2 CI band rather than at the top of the axes, where the
# band and its dashed line ran straight through the legend text.
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.012, ci_lo / ax.get_ylim()[1] - 0.02),
          fontsize=9.6, frameon=False, handlelength=2.3, labelspacing=1.0, borderpad=0.9)

table_lines = ['Chao2 at other clustering thresholds', '']
for row in rows:
    marker = ' ←' if row['threshold'] == PRIMARY else ''
    table_lines.append(f"{int(float(row['threshold']) * 100)}%   observed {row['observed']:>5,}   "
                       f"Chao2 {row['chao2']:>6,.0f}   {row['pct']:>2.0f}% sampled{marker}")
ax.text(0.985, 0.045, '\n'.join(table_lines), transform=ax.transAxes, ha='right', va='bottom',
        fontsize=8.6, color=TEXT_MUTED, family='DejaVu Sans Mono', linespacing=1.6)

fig.suptitle('PhaC cluster discovery has not saturated, and the rate depends on the study pool',
             fontsize=17, fontweight='bold', y=0.985, x=0.055, ha='left')
fig.text(0.055, 0.935,
         'Accumulation of distinct phaC sequence clusters, drawn two ways. Treating the study rather than the genome as the unit of independent sampling is the\n'
         'conservative estimate: the gap between the curves is diversity that pooling studies makes look cheaper to find. The two necessarily meet at the full sample,\n'
         'so the gap in between, not the endpoint, is the quantity of interest.',
         ha='left', va='top', fontsize=10.3, color=TEXT_MUTED, linespacing=1.5)
fig.tight_layout(rect=[0.005, 0.01, 0.995, 0.90])
fig.savefig(OUT / 'phac_accumulation_curves.png', dpi=300, facecolor='white')
fig.savefig(OUT / 'phac_accumulation_curves.pdf', facecolor='white')
print('saved', OUT / 'phac_accumulation_curves.png')
