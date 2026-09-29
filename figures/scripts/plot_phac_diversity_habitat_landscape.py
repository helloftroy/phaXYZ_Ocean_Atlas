"""'PhaC diversity accumulation landscape' -- replaces the four-panel
rank-abundance/rarefaction figure (plot_phac_rank_abundance_rarefaction.py)
with one hero figure that answers the same underlying question (is phaC
cluster diversity saturating, and does habitat explain any of it) more
directly: a single global rarefaction curve as the sampling-effort
baseline, with every well-sampled habitat plotted as one point against
that baseline, colored by whether it disproportionately contributes
globally rare clusters. Spec from the user (drafted with GPT, refined
here), reproduced in intent:

  1. Global rarefaction: repeatedly subsample N phaC-positive genomes at
     increasing N, count distinct phaC_cluster0.7 clusters, many
     permutations -> mean curve + 95% permutation interval. This is the
     richness EXPECTED from sampling effort alone, independent of
     habitat.
  2. Every habitat with >=50 phaC-positive genomes becomes one point:
     x = its genome count, y = its OBSERVED distinct-cluster count.
     Above the global curve = richer than expected for its sample size;
     below = more redundant.
  3. Point color = log2(observed rare-cluster fraction / expected
     rare-cluster fraction from matched random genome samples of the
     same size) -- "rare" defined GLOBALLY (a cluster found in <=5
     genomes anywhere in the atlas), not per-habitat, so a habitat can't
     be "rare" just by being small.
  4. Point size kept uniform throughout, per the user's own preference
     (avoids a busy/over-encoded figure).

Two technical points from the user, both load-bearing:
  - Richness is counted in DISTINCT CLUSTER IDs among the sampled/
    habitat genomes, never raw phaC protein/target_id count -- a genome
    with 2 paralogs in the same cluster contributes 1, matching
    plot_phac_rank_abundance_rarefaction.py's own convention exactly
    (same genome_clusters construction, reused verbatim below).
  - A habitat's rare-cluster fraction is computed over its DISTINCT
    clusters present, not occurrences -- a common cluster hit by 500
    genomes in one habitat still counts once.

Extrapolation (iNEXT/Chao-style) was requested only if straightforward,
with explicit instruction not to invent a curve shape if it isn't. A full
extrapolated RICHNESS CURVE (Colwell et al.'s unified rarefaction/
extrapolation framework) needs unsampled-frequency-count machinery this
project has no prior implementation of and no dependency (no iNEXT
equivalent in Python installed here) -- not attempted, per the
instruction. A single ASYMPTOTIC richness estimate (Chao2, incidence-
based -- Chao 1987) is a well-established, simple, one-line closed-form
statistic (not a curve), so that IS computed and reported as a headline
number alongside the observed rarefaction curve's own terminal marginal
discovery rate (slope over the last 5% of sampled genomes).

QC/exclusion conventions all match plot_phac_rank_abundance_rarefaction.py
exactly (same bad-target filter, same target_id->genome join fixing the
one-to-many bug documented there, same ecosystem_compartment habitat
field + EXCLUDE_HABITATS list) -- this figure is a replacement for that
one's panels A/C, not an independent re-derivation, so the underlying
population is deliberately identical.

Usage:
    python figures/scripts/plot_phac_diversity_habitat_landscape.py
"""
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

RANDOM_SEED = 42
N_GLOBAL_PERM = 300
N_HABITAT_PERM = 300
N_CURVE_POINTS = 200
MIN_HABITAT_N = 50
RARE_MAX_GENOMES = 5     # a cluster found in <=5 genomes globally is "rare" (singleton + rare tiers)
INTERMEDIATE_MAX_GENOMES = 20

EXCLUDE_HABITATS = {  # keep in sync with plot_phac_rank_abundance_rarefaction.py
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}

random.seed(RANDOM_SEED)
rng = np.random.default_rng(RANDOM_SEED)

# ---------------------------------------------------------------------
# 1. genome -> set(cluster_id), cluster_id -> set(genome)  (identical
#    construction to plot_phac_rank_abundance_rarefaction.py -- see that
#    script's own comments for why the one-to-many target_id->genome
#    join and bad-target filtering both matter here)
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()

target_to_genomes = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        if target_id in bad_targets:
            continue
        target_to_genomes[target_id].add(genome)

cluster_genomes = defaultdict(set)
genome_clusters = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        if target_id in bad_targets:
            continue
        for g in target_to_genomes.get(target_id, ()):
            cluster_genomes[cluster_id].add(g)
            genome_clusters[g].add(cluster_id)

all_genomes = np.array(sorted(genome_clusters))
n_genomes_total = len(all_genomes)
n_clusters_total = len(cluster_genomes)
print(f'{n_genomes_total:,} phaC-positive genomes, {n_clusters_total:,} distinct phaC_cluster0.7 clusters')

# global cluster-size tiers (by DISTINCT GENOMES carrying the cluster, project-wide)
cluster_size = {c: len(gs) for c, gs in cluster_genomes.items()}
n_singleton_global = sum(1 for s in cluster_size.values() if s == 1)
n_rare_global = sum(1 for s in cluster_size.values() if 2 <= s <= RARE_MAX_GENOMES)
n_intermediate_global = sum(1 for s in cluster_size.values() if RARE_MAX_GENOMES < s <= INTERMEDIATE_MAX_GENOMES)
n_common_global = sum(1 for s in cluster_size.values() if s > INTERMEDIATE_MAX_GENOMES)
rare_cluster_ids = {c for c, s in cluster_size.items() if s <= RARE_MAX_GENOMES}
print(f'global tiers: singleton={n_singleton_global:,} ({100*n_singleton_global/n_clusters_total:.1f}%), '
      f'rare(2-{RARE_MAX_GENOMES})={n_rare_global:,}, '
      f'intermediate({RARE_MAX_GENOMES+1}-{INTERMEDIATE_MAX_GENOMES})={n_intermediate_global:,}, '
      f'common(>{INTERMEDIATE_MAX_GENOMES})={n_common_global:,}')

# genome -> list of its cluster_ids, as a numpy-friendly structure for fast permutation
genome_cluster_list = [genome_clusters[g] for g in all_genomes]

# ---------------------------------------------------------------------
# 2. habitat assignment (same field/exclusion list as the figure this replaces)
# ---------------------------------------------------------------------
genome_habitat = {}
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c in EXCLUDE_HABITATS:
            continue
        genome_habitat[row['genome']] = c

habitat_genome_idx = defaultdict(list)
for i, g in enumerate(all_genomes):
    h = genome_habitat.get(g)
    if h:
        habitat_genome_idx[h].append(i)

valid_habitats = {h: idx for h, idx in habitat_genome_idx.items() if len(idx) >= MIN_HABITAT_N}
print(f'{len(valid_habitats)} habitats with >={MIN_HABITAT_N} phaC-positive genomes')

# ---------------------------------------------------------------------
# 3. global rarefaction: permutation mean + 95% percentile band
# ---------------------------------------------------------------------
def downsample_indices(n, n_points):
    if n <= n_points:
        return np.arange(1, n + 1)
    return np.unique(np.round(np.linspace(1, n, n_points)).astype(int))


xs_global = downsample_indices(n_genomes_total, N_CURVE_POINTS)
curves = np.zeros((N_GLOBAL_PERM, len(xs_global)), dtype=float)
order_template = np.arange(n_genomes_total)
for p in range(N_GLOBAL_PERM):
    order = rng.permutation(order_template)
    seen = set()
    xi = 0
    for i, gi in enumerate(order, start=1):
        seen.update(genome_cluster_list[gi])
        if xi < len(xs_global) and i == xs_global[xi]:
            curves[p, xi] = len(seen)
            xi += 1
print(f'global rarefaction: {N_GLOBAL_PERM} permutations done')

mean_curve = curves.mean(axis=0)
lo_curve = np.percentile(curves, 2.5, axis=0)
hi_curve = np.percentile(curves, 97.5, axis=0)

# terminal marginal discovery rate: slope of the mean curve over its last 5%
tail_start = int(0.95 * len(xs_global))
dx = xs_global[-1] - xs_global[tail_start]
dy = mean_curve[-1] - mean_curve[tail_start]
marginal_rate_per_1000 = 1000 * dy / dx if dx > 0 else float('nan')
print(f'terminal marginal discovery rate: {marginal_rate_per_1000:.2f} new clusters per 1,000 additional genomes '
      f'(measured over the last {dx:,} of {n_genomes_total:,} genomes)')

# Chao2 asymptotic richness estimate (incidence-based, Chao 1987) -- a
# single number, not a curve. Q1 = clusters seen in exactly 1 genome
# (= n_singleton_global), Q2 = clusters seen in exactly 2 genomes.
Q1 = n_singleton_global
Q2 = sum(1 for s in cluster_size.values() if s == 2)
m = n_genomes_total
if Q2 > 0:
    chao2 = n_clusters_total + ((m - 1) / m) * (Q1 ** 2) / (2 * Q2)
else:
    chao2 = n_clusters_total + ((m - 1) / m) * (Q1 * (Q1 - 1)) / 2
print(f'Chao2 asymptotic richness estimate: {chao2:,.0f} (vs. {n_clusters_total:,} observed, '
      f'{100 * n_clusters_total / chao2:.1f}% of estimated total already sampled)')

# ---------------------------------------------------------------------
# 4. per-habitat: observed richness/rare-fraction + matched-n null permutation
# ---------------------------------------------------------------------
habitat_stats = []
for h, idx in sorted(valid_habitats.items(), key=lambda kv: -len(kv[1])):
    idx = np.array(idx)
    n_h = len(idx)
    observed_set = set()
    for i in idx:
        observed_set |= genome_cluster_list[i]
    n_obs = len(observed_set)
    n_singleton_h = sum(1 for c in observed_set if cluster_size[c] == 1)
    n_rare_h = sum(1 for c in observed_set if 2 <= cluster_size[c] <= RARE_MAX_GENOMES)
    rare_fraction_obs = (n_singleton_h + n_rare_h) / n_obs if n_obs else 0.0

    null_richness = np.zeros(N_HABITAT_PERM)
    null_rare_frac = np.zeros(N_HABITAT_PERM)
    for p in range(N_HABITAT_PERM):
        draw = rng.choice(order_template, size=n_h, replace=False)
        draw_set = set()
        for gi in draw:
            draw_set |= genome_cluster_list[gi]
        null_richness[p] = len(draw_set)
        n_rare_draw = sum(1 for c in draw_set if cluster_size[c] <= RARE_MAX_GENOMES)
        null_rare_frac[p] = n_rare_draw / len(draw_set) if draw_set else 0.0

    expected = null_richness.mean()
    ci_lo, ci_hi = np.percentile(null_richness, [2.5, 97.5])
    richness_ratio = n_obs / expected if expected else float('nan')
    richness_residual = n_obs - expected
    p_richness = ((null_richness >= n_obs).mean() if n_obs >= expected
                  else (null_richness <= n_obs).mean())

    expected_rare_frac = null_rare_frac.mean()
    rare_ci_lo, rare_ci_hi = np.percentile(null_rare_frac, [2.5, 97.5])
    rare_enrichment = (np.log2(rare_fraction_obs / expected_rare_frac)
                       if rare_fraction_obs > 0 and expected_rare_frac > 0
                       else (float('-inf') if expected_rare_frac > 0 else float('nan')))
    p_rare = ((null_rare_frac >= rare_fraction_obs).mean() if rare_fraction_obs >= expected_rare_frac
              else (null_rare_frac <= rare_fraction_obs).mean())

    habitat_stats.append(dict(
        habitat=h, n_genomes=n_h, n_clusters=n_obs,
        expected_clusters=expected, expected_clusters_ci_lo=ci_lo, expected_clusters_ci_hi=ci_hi,
        richness_ratio=richness_ratio, richness_residual=richness_residual, richness_p=p_richness,
        n_singleton_clusters=n_singleton_h, n_rare_clusters=n_rare_h,
        rare_fraction=rare_fraction_obs, expected_rare_fraction=expected_rare_frac,
        expected_rare_fraction_ci_lo=rare_ci_lo, expected_rare_fraction_ci_hi=rare_ci_hi,
        rare_enrichment=rare_enrichment, rare_p=p_rare,
    ))
    print(f'  {h:32s} n={n_h:>6,} obs_clusters={n_obs:>5,} expected={expected:>7.1f} '
          f'ratio={richness_ratio:.2f} rare_enrich={rare_enrichment:+.2f}')

# ---------------------------------------------------------------------
# 5. write stats table
# ---------------------------------------------------------------------
stats_path = OUT / 'phac_diversity_habitat_landscape_stats.tsv'
fieldnames = ['habitat', 'n_genomes', 'n_clusters', 'expected_clusters', 'expected_clusters_ci_lo',
              'expected_clusters_ci_hi', 'richness_ratio', 'richness_residual', 'richness_p',
              'n_singleton_clusters', 'n_rare_clusters', 'rare_fraction', 'expected_rare_fraction',
              'expected_rare_fraction_ci_lo', 'expected_rare_fraction_ci_hi', 'rare_enrichment', 'rare_p']
with open(stats_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    for row in habitat_stats:
        writer.writerow(row)
print(f'wrote {stats_path}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5})

# diverging colormap, validated CVD-safe pair (dataviz skill palette):
# blue #2a78d6 (cool, negative/depleted) <-> gray #f0efec (neutral) <-> red #e34948 (warm, enriched)
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    'rare_enrichment', ['#2a78d6', '#f0efec', '#e34948'], N=256)
CURVE_COLOR = '#3A4A46'      # neutral charcoal -- the baseline, not a competing hue
CURVE_BAND = '#E4E8E5'
TEXT_MUTED = '#5B6E70'
EDGE_WHITE = '#FCFCFB'

fig, ax = plt.subplots(figsize=(14, 9.8), dpi=300)

# global expectation curve + 95% permutation band
ax.fill_between(xs_global, lo_curve, hi_curve, color=CURVE_BAND, zorder=1, linewidth=0,
                 label='95% permutation interval')
ax.plot(xs_global, mean_curve, color=CURVE_COLOR, linewidth=2.2, zorder=2,
        label='Global expectation (mean of\nrandom genome subsamples)')

# habitat points
color_vals = np.array([hs['rare_enrichment'] for hs in habitat_stats])
finite_vals = color_vals[np.isfinite(color_vals)]
vlim = np.nanpercentile(np.abs(finite_vals), 95) if len(finite_vals) else 1.0
vlim = max(vlim, 0.5)
clipped_vals = np.clip(np.nan_to_num(color_vals, nan=0.0, neginf=-vlim, posinf=vlim), -vlim, vlim)

xs_pts = np.array([hs['n_genomes'] for hs in habitat_stats])
ys_pts = np.array([hs['n_clusters'] for hs in habitat_stats])
sc = ax.scatter(xs_pts, ys_pts, c=clipped_vals, cmap=DIVERGING_CMAP, vmin=-vlim, vmax=vlim,
                 s=210, edgecolors=EDGE_WHITE, linewidths=1.6, zorder=5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('PhaC-positive genomes sampled (log scale)', fontsize=12)
ax.set_ylabel('Distinct phaC_cluster0.7 clusters recovered (log scale)', fontsize=12)
ax.grid(True, which='major', color='#EBEEEC', linewidth=0.6, zorder=0)
ax.grid(True, which='minor', color='#F4F6F4', linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_color('#C3C2B7')
ax.spines['bottom'].set_color('#C3C2B7')

# selective labels. Direct point labels only for habitats with enough
# surrounding white space to read cleanly -- most of this dataset's
# n>=50 habitats sit in one tight n=100-300 clump (see the crowded-region
# crop checked during design), where any inline label collides with
# neighboring points. Those instead go in a side callout list (below),
# linked to their point only by matching color, not a leader line into
# the clutter -- direct labeling is preferred throughout this project's
# figures, but only where it stays legible; forcing it into a dense
# cluster does the opposite of the "selective" instruction.
DENSE_ZONE_MAX_N = 400
by_ratio = sorted(habitat_stats, key=lambda h: h['richness_ratio'])
by_rare = sorted([h for h in habitat_stats if np.isfinite(h['rare_enrichment'])], key=lambda h: h['rare_enrichment'])
largest_n = max(habitat_stats, key=lambda h: h['n_genomes'])

direct_label_candidates = {h['habitat']: h for h in (by_ratio[:1] + by_ratio[-1:] + [largest_n])}
direct_labels = {name: h for name, h in direct_label_candidates.items() if h['n_genomes'] > DENSE_ZONE_MAX_N}

label_offsets = {
    'Seawater': (16, 8),
    'Marine Porifera tissue': (18, -6),
    'Marine sediment': (-20, 26),
}
for name, h in direct_labels.items():
    dx_lab, dy_lab = label_offsets.get(name, (16, 10))
    ax.annotate(
        name, xy=(h['n_genomes'], h['n_clusters']), xytext=(dx_lab, dy_lab),
        textcoords='offset points', fontsize=10.5, color='#2A2E2C', fontweight='medium',
        path_effects=[pe.withStroke(linewidth=3, foreground=EDGE_WHITE)],
        arrowprops=dict(arrowstyle='-', color='#9A9890', linewidth=0.9, shrinkA=7, shrinkB=7), zorder=6)

# side callout: top 3 rare-cluster enrichments, almost all of which fall
# inside the dense n<400 clump -- named here instead of inline.
top_rare = [h for h in reversed(by_rare) if h['habitat'] not in direct_labels][:3]
callout_lines = [(h['habitat'], h['rare_enrichment'], DIVERGING_CMAP((h['rare_enrichment'] + vlim) / (2 * vlim)))
                 for h in top_rare]

# render the callout as small color swatches + text, stacked in the
# open lower-right area under the curve (no leader lines needed).
callout_x0, callout_y0 = 0.615, 0.30
ax.text(callout_x0, callout_y0, 'Highest rare-cluster enrichment', transform=ax.transAxes,
        ha='left', va='top', fontsize=10, fontweight='bold', color='#2A2E2C')
for i, (name, val, color) in enumerate(callout_lines):
    y = callout_y0 - 0.045 * (i + 1.15)
    ax.scatter([callout_x0 + 0.012], [y], transform=ax.transAxes, s=130, color=color,
               edgecolors=EDGE_WHITE, linewidths=1.4, clip_on=False, zorder=7)
    ax.text(callout_x0 + 0.032, y, f'{name}  (+{val:.1f})', transform=ax.transAxes,
            ha='left', va='center', fontsize=9.7, color='#2A2E2C')

# colorbar (diverging, centered at 0 -- the "legend" for the continuous encoding)
cbar = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045, aspect=26)
cbar.set_label('Rare-cluster enrichment  log$_2$(observed / expected)', fontsize=10.5, color='#3A3E3B')
cbar.ax.tick_params(labelsize=9.5)
cbar.outline.set_visible(False)

ax.legend(loc='lower right', fontsize=9.5, frameon=False, handlelength=2.2)

# headline stats box
stats_text = (
    f'{n_clusters_total:,} clusters from {n_genomes_total:,} phaC-positive genomes\n'
    f'{n_singleton_global:,} singletons ({100*n_singleton_global/n_clusters_total:.0f}%) -- '
    f'rare ($\\leq${RARE_MAX_GENOMES} genomes): {n_singleton_global+n_rare_global:,} ({100*(n_singleton_global+n_rare_global)/n_clusters_total:.0f}%)\n'
    f'Chao2 estimated total richness: {chao2:,.0f}  ({100*n_clusters_total/chao2:.0f}% already observed)\n'
    f'Terminal discovery rate: {marginal_rate_per_1000:.1f} new clusters / 1,000 genomes'
)
ax.text(0.015, 0.985, stats_text, transform=ax.transAxes, ha='left', va='top', fontsize=9.3,
        color='#3A4A46', linespacing=1.55,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F7F8F6', edgecolor='#D8DCD9', linewidth=0.8))

fig.suptitle('The PhaC diversity accumulation landscape', fontsize=19, fontweight='bold', y=0.985, x=0.065, ha='left')
fig.text(0.065, 0.945,
         'Habitats above the curve hold more phaC diversity than sampling effort alone predicts; below, less.\n'
         'Color shows enrichment for globally rare clusters, defined independent of each habitat’s own size.',
         ha='left', va='top', fontsize=10.5, color=TEXT_MUTED, linespacing=1.5)

fig.tight_layout(rect=[0.005, 0.01, 0.99, 0.90])

out_path = OUT / 'phac_diversity_habitat_landscape.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_diversity_habitat_landscape.pdf', facecolor='white')
