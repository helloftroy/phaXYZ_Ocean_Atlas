"""Community-ecology framing for the phaC_cluster0.7 cluster set: rank-
abundance (a few cosmopolitan clusters, long tail of rare/singleton ones)
and discovery/rarefaction curves (does cumulative sampling still turn up
new clusters, or has it flattened out). The point, per the user's own
framing, is honesty about scope: this dataset is a first atlas of phaC
diversity, built from whatever metagenomes happen to be public, not a
saturated/complete survey -- these curves are the direct evidence for
that claim (or against it, if they show early saturation).

Cluster "abundance" = DISTINCT GENOMES carrying that cluster, not raw
target_id (protein hit) count -- a genome with two phaC paralogs in the
same cluster should count once, the same way a community ecologist counts
individuals, not gene copies. Built directly from the raw mmseqs2 cluster
membership file (phaC_cluster0.7_cluster.tsv, 20,211 clusters) joined to
the clean target_id->genome table (phaC_all_genomes_from_nr100_clusters.tsv),
rather than reusing the pre-aggregated phaC_cluster0.7_cluster_ecology.tsv
(which only covers 11,776/20,211 clusters after an unclear contamination-
filtering pass -- computing this directly from the raw membership file
keeps every cluster in scope and the provenance auditable).

Four panels:
  A. Rank-abundance (log-log): all 20,211 clusters, ranked by n_genomes.
  B. Discovery curve, genome-order vs. study-order shuffling: random
     genome order will look closer to "saturated" than reality, because a
     handful of very large individual studies (top 5 alone = 36% of all
     phaC-positive genomes -- confirmed directly below) pad the curve with
     many genomes from places already sampled. Study-block shuffling
     (whole studies added at once, in random study order) removes that
     bias -- if it still climbs near the end, that's a fairer signal that
     new studies keep finding new clusters, not just resampling known ones.
  C. Discovery curve, well-sampled vs. rare/host-associated habitats: same
     ≥200-genome habitat set as PHA_CLEAN_RESULTS.md §5.3, split at a
     natural order-of-magnitude gap in total screened genomes (>10,000 vs
     <5,000) into "well-sampled" (Seawater, Marine sediment, Sponge tissue)
     and "rare/host-associated/extreme" (the other 12: whale-fall bone
     biofilm, hydrothermal vent categories, cold seep sediment, sea ice,
     coral/hydrozoa/algae tissue, estuarine categories, brackish water,
     biofilm). Plotted on a shared log-x (raw genome count) axis so the
     SLOPE at each curve's own endpoint is the thing to read: a curve
     still rising steeply at its right edge means sampling more genomes
     from that bucket would likely keep finding new clusters; a flat
     right edge means that bucket is closer to saturated.
  D. Same idea, shallow vs. deep, using real (NCBI BioSample-derived)
     depth_zone -- caveated hard, because coverage is genuinely sparse
     (see below) and there is no dedicated "hadal" bucket in the data at
     all, only a >4000m catch-all with a small n. Reported honestly as a
     secondary, lower-confidence check, not the headline panel.

Depth coverage caveat (checked directly, not assumed): genome_family_matrix.tsv's
own depth_raw/depth_m/depth_zone columns are entirely empty for every row --
the real depth data lives only in phaC_unique_targets_with_metadata_depth.tsv
(NCBI BioSample-derived, see phaatlas/pipeline/ncbi_depth.py), and even
there only ~33% of target rows have any depth_zone at all, with the
deepest bucket (>4000m) covering roughly 1% of all rows. Panel D is
built from that sparse subset only and says so on the figure itself.

Usage:
    python figures/scripts/plot_phac_rank_abundance_rarefaction.py
"""
import csv
import random
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

RANDOM_SEED = 42
N_PERMUTATIONS = 30       # for genome-order and habitat/depth-bucket curves
N_STUDY_PERMUTATIONS = 200  # study-block shuffles are cheap (only ~195 blocks) so run more
N_CURVE_POINTS = 250       # downsample each curve to this many points for plotting/averaging

EXCLUDE_HABITATS = {  # keep in sync with plot_phac_pct_by_ocean_habitat.py
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
MIN_HABITAT_N = 200
WELL_SAMPLED_HABITATS = {'Seawater', 'Marine sediment', 'Marine Porifera tissue'}  # screened n > 10,000

random.seed(RANDOM_SEED)

# ---------------------------------------------------------------------
# 1. cluster_id -> set(genome), genome -> set(cluster_id)   [full 20,211 clusters]
# ---------------------------------------------------------------------
target_to_genome = {}
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)  # header
    for target_id, genome in r:
        target_to_genome[target_id] = genome

cluster_genomes = defaultdict(set)
genome_clusters = defaultdict(set)
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        g = target_to_genome.get(target_id)
        if g is None:
            continue
        cluster_genomes[cluster_id].add(g)
        genome_clusters[g].add(cluster_id)

n_clusters = len(cluster_genomes)
all_genomes = sorted(genome_clusters)  # sorted first for determinism, shuffled later
n_genomes_total = len(all_genomes)
print(f'{n_clusters:,} clusters, {n_genomes_total:,} distinct genomes (from {len(target_to_genome):,} target_id->genome rows)')

sizes = sorted((len(gs) for gs in cluster_genomes.values()), reverse=True)
n_singletons = sum(1 for s in sizes if s == 1)
print(f'largest cluster: {sizes[0]:,} genomes; {n_singletons:,}/{n_clusters:,} clusters ({100*n_singletons/n_clusters:.1f}%) are singletons (1 genome)')

# ---------------------------------------------------------------------
# 2. genome -> study_id  (for study-block shuffling)
# ---------------------------------------------------------------------
genome_study = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['study_id']:
            genome_study[row['genome']] = row['study_id']

genomes_with_study = [g for g in all_genomes if g in genome_study]
print(f'{len(genomes_with_study):,}/{n_genomes_total:,} phaC-cluster genomes have a study_id')
study_genomes = defaultdict(list)
for g in genomes_with_study:
    study_genomes[genome_study[g]].append(g)
top5 = sorted(study_genomes.values(), key=len, reverse=True)[:5]
top5_share = 100 * sum(len(v) for v in top5) / len(genomes_with_study)
print(f'{len(study_genomes)} studies; top 5 alone = {top5_share:.1f}% of genomes with a study_id')

# ---------------------------------------------------------------------
# 3. genome -> ecosystem_compartment (fresh canonical habitat labels, same
#    file/convention as plot_phac_pct_by_ocean_habitat.py)
# ---------------------------------------------------------------------
genome_habitat = {}
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c in EXCLUDE_HABITATS:
            continue
        genome_habitat[row['genome']] = c

habitat_totals = defaultdict(int)
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c not in EXCLUDE_HABITATS:
            habitat_totals[c] += 1
valid_habitats = {c for c, n in habitat_totals.items() if n >= MIN_HABITAT_N}

well_sampled_genomes = [g for g in all_genomes if genome_habitat.get(g) in WELL_SAMPLED_HABITATS]
rare_habitat_genomes = [g for g in all_genomes
                         if genome_habitat.get(g) in (valid_habitats - WELL_SAMPLED_HABITATS)]
print(f'well-sampled habitat genomes (phaC-positive): {len(well_sampled_genomes):,}')
print(f'rare/host-associated/extreme habitat genomes (phaC-positive): {len(rare_habitat_genomes):,}')

# ---------------------------------------------------------------------
# 4. genome -> depth_zone (sparse, NCBI-BioSample-derived; see module docstring)
# ---------------------------------------------------------------------
genome_depth_zone = {}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['depth_zone'] and row['genome'] not in genome_depth_zone:
            genome_depth_zone[row['genome']] = row['depth_zone']

SHALLOW_ZONES = {'0-50m', '50-200m'}
DEEP_ZONES = {'200-1000m', '1000-4000m', '>4000m'}
shallow_genomes = [g for g in all_genomes if genome_depth_zone.get(g) in SHALLOW_ZONES]
deep_genomes = [g for g in all_genomes if genome_depth_zone.get(g) in DEEP_ZONES]
n_with_depth = sum(1 for g in all_genomes if g in genome_depth_zone)
print(f'{n_with_depth:,}/{n_genomes_total:,} genomes ({100*n_with_depth/n_genomes_total:.1f}%) have any depth_zone; '
      f'shallow={len(shallow_genomes):,}, deep={len(deep_genomes):,}')


# ---------------------------------------------------------------------
# discovery-curve machinery
# ---------------------------------------------------------------------
def downsample_indices(n, n_points):
    if n <= n_points:
        return list(range(1, n + 1))
    idx = np.unique(np.round(np.linspace(1, n, n_points)).astype(int))
    return idx.tolist()


def discovery_curve_genome_order(genome_list, n_perm):
    """Mean cumulative distinct-cluster count at each sampled genome-count,
    averaged over n_perm random shuffles of genome_list."""
    n = len(genome_list)
    if n == 0:
        return np.array([]), np.array([]), np.array([])
    xs = downsample_indices(n, N_CURVE_POINTS)
    all_curves = np.zeros((n_perm, len(xs)), dtype=float)
    for p in range(n_perm):
        order = genome_list[:]
        random.shuffle(order)
        seen = set()
        curve = []
        xi = 0
        for i, g in enumerate(order, start=1):
            seen.update(genome_clusters[g])
            if xi < len(xs) and i == xs[xi]:
                curve.append(len(seen))
                xi += 1
        all_curves[p, :len(curve)] = curve
    mean = all_curves.mean(axis=0)
    lo = all_curves.min(axis=0)
    hi = all_curves.max(axis=0)
    return np.array(xs), mean, (lo, hi)


def discovery_curve_study_order(study_genomes_map, n_perm):
    """Same idea but shuffles at the STUDY level -- each study's full genome
    list is added as one block, in random study order, so within-study
    genome order never matters (removes the 'a few huge studies pad the
    curve' bias flagged in the module docstring)."""
    studies = list(study_genomes_map.keys())
    total_n = sum(len(v) for v in study_genomes_map.values())
    xs = downsample_indices(total_n, N_CURVE_POINTS)
    all_curves = np.zeros((n_perm, len(xs)), dtype=float)
    for p in range(n_perm):
        order = studies[:]
        random.shuffle(order)
        seen = set()
        curve = []
        xi = 0
        i = 0
        for s in order:
            for g in study_genomes_map[s]:
                i += 1
                seen.update(genome_clusters[g])
                if xi < len(xs) and i == xs[xi]:
                    curve.append(len(seen))
                    xi += 1
        while xi < len(xs):  # final block may overshoot past the last few downsample points
            curve.append(len(seen))
            xi += 1
        all_curves[p, :len(curve)] = curve[:len(xs)]
    mean = all_curves.mean(axis=0)
    lo = all_curves.min(axis=0)
    hi = all_curves.max(axis=0)
    return np.array(xs), mean, (lo, hi)


# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
axA, axB, axC, axD = axes.flat

# Panel A: rank-abundance
ranks = np.arange(1, len(sizes) + 1)
axA.plot(ranks, sizes, color='#1E6E7A', linewidth=1.6)
axA.set_xscale('log')
axA.set_yscale('log')
axA.set_xlabel('Cluster rank (log)')
axA.set_ylabel('Genomes in cluster (log)')
axA.set_title('A. Rank-abundance of phaC clusters', fontsize=12.5, fontweight='bold', loc='left')
axA.text(0.97, 0.95, f'{n_clusters:,} clusters\n{n_singletons:,} singletons ({100*n_singletons/n_clusters:.0f}%)\nlargest = {sizes[0]:,} genomes',
          transform=axA.transAxes, ha='right', va='top', fontsize=9, color='#3A4A46',
          bbox=dict(boxstyle='round', facecolor='#F2F5F3', edgecolor='#C8D2CD'))
axA.grid(True, which='both', color='#E4E8E5', linewidth=0.5, zorder=0)
axA.set_axisbelow(True)

# Panel B: genome-order vs study-order discovery
xs_g, mean_g, (lo_g, hi_g) = discovery_curve_genome_order(genomes_with_study, N_PERMUTATIONS)
xs_s, mean_s, (lo_s, hi_s) = discovery_curve_study_order(study_genomes, N_STUDY_PERMUTATIONS)
axB.fill_between(xs_g, lo_g, hi_g, color='#1E6E7A', alpha=0.18, linewidth=0)
axB.plot(xs_g, mean_g, color='#1E6E7A', linewidth=1.8, label='Random genome order')
axB.fill_between(xs_s, lo_s, hi_s, color='#C2A83E', alpha=0.22, linewidth=0)
axB.plot(xs_s, mean_s, color='#C2A83E', linewidth=1.8, label='Random study-block order')
axB.set_xlabel('Cumulative genomes sampled')
axB.set_ylabel('Distinct clusters discovered')
axB.set_title('B. Discovery curve: genome order vs. study order', fontsize=12.5, fontweight='bold', loc='left')
axB.legend(loc='lower right', fontsize=9, frameon=False)
axB.text(0.03, 0.95, f'top 5 studies = {top5_share:.0f}% of genomes\n({len(study_genomes)} studies total)',
          transform=axB.transAxes, ha='left', va='top', fontsize=8.7, color='#5B6E70')
axB.grid(True, color='#E4E8E5', linewidth=0.6, zorder=0)
axB.set_axisbelow(True)

# Panel C: well-sampled vs rare/host-associated habitats
xs_w, mean_w, (lo_w, hi_w) = discovery_curve_genome_order(well_sampled_genomes, N_PERMUTATIONS)
xs_r, mean_r, (lo_r, hi_r) = discovery_curve_genome_order(rare_habitat_genomes, N_PERMUTATIONS)
axC.fill_between(xs_w, lo_w, hi_w, color='#1E6E7A', alpha=0.18, linewidth=0)
axC.plot(xs_w, mean_w, color='#1E6E7A', linewidth=1.8,
          label=f'Well-sampled (Seawater, sediment, sponge; n={len(well_sampled_genomes):,})')
axC.fill_between(xs_r, lo_r, hi_r, color='#9E3B3B', alpha=0.18, linewidth=0)
axC.plot(xs_r, mean_r, color='#9E3B3B', linewidth=1.8,
          label=f'Rare/host-associated/extreme (12 habitats; n={len(rare_habitat_genomes):,})')
axC.set_xscale('log')
axC.set_xlabel('Cumulative genomes sampled (log)')
axC.set_ylabel('Distinct clusters discovered')
axC.set_title('C. Discovery curve by habitat sampling depth', fontsize=12.5, fontweight='bold', loc='left')
axC.legend(loc='lower right', fontsize=8.3, frameon=False)
axC.grid(True, which='both', color='#E4E8E5', linewidth=0.5, zorder=0)
axC.set_axisbelow(True)

# Panel D: shallow vs deep (sparse depth data -- caveated on the panel itself)
xs_sh, mean_sh, (lo_sh, hi_sh) = discovery_curve_genome_order(shallow_genomes, N_PERMUTATIONS)
xs_dp, mean_dp, (lo_dp, hi_dp) = discovery_curve_genome_order(deep_genomes, N_PERMUTATIONS)
axD.fill_between(xs_sh, lo_sh, hi_sh, color='#1E6E7A', alpha=0.18, linewidth=0)
axD.plot(xs_sh, mean_sh, color='#1E6E7A', linewidth=1.8, label=f'Shallow, 0-200m (n={len(shallow_genomes):,})')
axD.fill_between(xs_dp, lo_dp, hi_dp, color='#9E3B3B', alpha=0.18, linewidth=0)
axD.plot(xs_dp, mean_dp, color='#9E3B3B', linewidth=1.8, label=f'Deep, 200->4000m (n={len(deep_genomes):,})')
axD.set_xlabel('Cumulative genomes sampled')
axD.set_ylabel('Distinct clusters discovered')
axD.set_title('D. Discovery curve by depth (sparse data -- see caveat)', fontsize=12.5, fontweight='bold', loc='left')
axD.legend(loc='lower right', fontsize=9, frameon=False)
axD.text(0.03, 0.95, f'only {100*n_with_depth/n_genomes_total:.0f}% of genomes have any depth_zone;\nno dedicated hadal bucket exists in this data',
          transform=axD.transAxes, ha='left', va='top', fontsize=8, color='#9E3B3B', style='italic')
axD.grid(True, color='#E4E8E5', linewidth=0.6, zorder=0)
axD.set_axisbelow(True)

for ax in axes.flat:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle('phaC diversity: rank-abundance and sampling-effort framing', fontsize=16, fontweight='bold', y=0.995)
fig.text(0.5, 0.965,
          f'{n_clusters:,} phaC_cluster0.7 clusters (70% identity) from {n_genomes_total:,} phaC-positive genomes. '
          f'Shaded bands = min-max across {N_PERMUTATIONS} (genome-order) / {N_STUDY_PERMUTATIONS} (study-order) random-order permutations.',
          ha='center', fontsize=9.3, color='#5B6E70')
fig.tight_layout(rect=[0, 0, 1, 0.955])

out_path = OUT / 'phac_rank_abundance_rarefaction.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_rank_abundance_rarefaction.pdf', facecolor='white')
