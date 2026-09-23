"""Does phaC sequence identity (phaC_cluster0.7 cluster membership) pair up
with which mcl-PHA precursor route a genome uses -- phaJ (beta-oxidation-
linked) vs. phaG (FAS-linked) vs. both vs. neither? The biological
motivation: a synthase that will receive its monomer from a fatty-acid-
synthesis-linked route vs. a beta-oxidation-linked route is plausibly under
different selective pressure (different upstream partner enzymes, possibly
different substrate-channeling requirements), so a real pairing would show
up as specific phaC clusters being disproportionately associated with one
route rather than the two routes being scattered evenly across every
cluster regardless of which phaC variant a genome carries.

RESTRICTED TO SINGLE-COPY GENOMES ONLY (n_phaC==1, from the same
authoritative genome_family_matrix.tsv column used in
plot_phac_multicopy_genomes.py). This is not an arbitrary simplification --
it is what makes the question well-posed at all: a genome with 2+ phaC
copies could pair ANY of its synthases with its one shared precursor-route
gene set, so "this cluster pairs with G" would be ambiguous by construction
for multi-copy genomes. Single-copy genomes give one phaC sequence, one
cluster, one unambiguous genome-level G/J read.

Two independent pipelines' single-copy calls are cross-checked rather than
trusted blindly: genome_family_matrix.tsv's n_phaC==1 (the family-level HMM/
cluster-derived gene count) vs. the raw NR100 cluster-membership join
(phaC_cluster0.7_cluster.tsv joined through phaC_all_genomes_from_nr100_clusters.tsv,
which -- unlike a naive one-to-one join -- is a genuine one-to-many mapping,
since its target_id is a 100%-identity-clustered REPRESENTATIVE shared by
every genome carrying an identical sequence). A genome only enters this
analysis if BOTH pipelines agree it has exactly one phaC target_id; genomes
where they disagree are excluded and counted, not silently resolved one way.

Usage:
    python figures/scripts/plot_phac_cluster_vs_precursor_route.py

Outputs:
    figures/phac_cluster_vs_precursor_route.png / .pdf
    figures/phac_cluster_precursor_route_all_clusters.tsv   (every qualifying cluster, sorted by phaG rate)
"""
import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats_utils import mantel_haenszel

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

MIN_CLUSTER_N = 30    # clusters smaller than this give an unstable %phaG estimate
TOP_N_CLUSTERS_PLOT = 15

# ---------------------------------------------------------------------
# 1. per-genome n_phaC, G/J presence, phylum (authoritative, phaC-positive only)
# ---------------------------------------------------------------------
genome_info = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genome_info[row['genome']] = {
            'n_phac': n,
            'g': int(row['n_phaG']) > 0,
            'j': int(row['n_phaJ']) > 0,
            'phylum': row['gtdb_phylum'],
        }

single_copy_a = {g for g, r in genome_info.items() if r['n_phac'] == 1}
print(f'{len(genome_info):,} phaC-positive genomes; {len(single_copy_a):,} single-copy per genome_family_matrix.tsv (n_phaC==1)')

# ---------------------------------------------------------------------
# 2. genome -> set(target_id) from the raw NR100 join (properly one-to-many
#    in BOTH directions -- see module docstring)
# ---------------------------------------------------------------------
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        genome_targets[genome].add(target_id)

single_copy_b = {g for g, targets in genome_targets.items() if len(targets) == 1}

single_copy_genomes = single_copy_a & single_copy_b
n_conflict = len(single_copy_a - single_copy_b) + len(single_copy_b - single_copy_a)
print(f'{len(single_copy_b):,} single-target per the raw NR100 join; '
      f'{len(single_copy_genomes):,} genomes agree by BOTH pipelines (used below); '
      f'{n_conflict:,} disagree or are missing from one pipeline (excluded, not resolved either way)')

# ---------------------------------------------------------------------
# 3. target_id -> cluster_id (confirmed 1:1, see script history) + genome -> cluster_id
#    for the clean single-copy set only
# ---------------------------------------------------------------------
target_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    for cluster_id, target_id in r:
        target_cluster[target_id] = cluster_id

genome_cluster = {}
for g in single_copy_genomes:
    (target_id,) = genome_targets[g]
    c = target_cluster.get(target_id)
    if c is not None:
        genome_cluster[g] = c

n_final = len(genome_cluster)
print(f'{n_final:,} single-copy genomes have a resolved phaC_cluster0.7 cluster -- this is the analysis set')


def mcl_strategy(rec):
    if rec['g'] and rec['j']:
        return 'Both (G+J)'
    if rec['g']:
        return 'FAS-linked only (G)'
    if rec['j']:
        return 'Beta-oxidation only (J)'
    return 'Neither'


baseline_n_g = sum(1 for g in genome_cluster if genome_info[g]['g'])
baseline_pct_g = 100 * baseline_n_g / n_final
print(f'baseline %phaG among these {n_final:,} single-copy genomes: {baseline_pct_g:.2f}%')

# ---------------------------------------------------------------------
# 4. per-cluster composition + phylum-controlled CMH test for phaG,
#    for every cluster with >= MIN_CLUSTER_N single-copy genomes
# ---------------------------------------------------------------------
cluster_genomes = defaultdict(list)
for g, c in genome_cluster.items():
    cluster_genomes[c].append(g)

qualifying_clusters = [c for c, gs in cluster_genomes.items() if len(gs) >= MIN_CLUSTER_N]
print(f'{len(qualifying_clusters)} clusters with >= {MIN_CLUSTER_N} single-copy genomes (of {len(cluster_genomes)} total touched)')

cluster_rows = []
for c in qualifying_clusters:
    gs = cluster_genomes[c]
    n = len(gs)
    recs = [genome_info[g] for g in gs]
    n_g_c = sum(1 for r in recs if r['g'])
    n_j_c = sum(1 for r in recs if r['j'])
    strat = Counter(mcl_strategy(r) for r in recs)
    top_phylum, top_phylum_n = Counter(r['phylum'] for r in recs).most_common(1)[0]

    # phylum-controlled CMH: is phaG still associated with THIS cluster (vs.
    # every other cluster in the single-copy analysis set) once phylum is
    # controlled -- same reasoning/machinery as the habitat tests.
    phylum_in = defaultdict(lambda: [0, 0])
    phylum_out = defaultdict(lambda: [0, 0])
    for g, c2 in genome_cluster.items():
        rec = genome_info[g]
        bucket = phylum_in if c2 == c else phylum_out
        bucket[rec['phylum']][0 if rec['g'] else 1] += 1
    strata = [(*phylum_in.get(p, [0, 0]), *phylum_out.get(p, [0, 0])) for p in set(phylum_in) | set(phylum_out)]
    or_mh, cmh_stat, cmh_p, n_strata = mantel_haenszel(strata)

    cluster_rows.append({
        'cluster_id': c, 'n_genomes': n, 'pct_phaG': f'{100*n_g_c/n:.2f}', 'pct_phaJ': f'{100*n_j_c/n:.2f}',
        'n_neither': strat.get('Neither', 0), 'n_beta_oxidation_only_J': strat.get('Beta-oxidation only (J)', 0),
        'n_fas_linked_only_G': strat.get('FAS-linked only (G)', 0), 'n_both_GJ': strat.get('Both (G+J)', 0),
        'top_phylum': top_phylum, 'top_phylum_pct': f'{100*top_phylum_n/n:.1f}',
        'cmh_or_phaG': f'{or_mh:.3f}' if or_mh is not None else '',
        'cmh_p_phaG': f'{cmh_p:.3g}' if cmh_p is not None else '',
        'n_strata_used': n_strata,
    })

cluster_rows.sort(key=lambda r: -float(r['pct_phaG']))
cluster_out = OUT / 'phac_cluster_precursor_route_all_clusters.tsv'
with open(cluster_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(cluster_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(cluster_rows)
print(f'wrote {cluster_out} ({len(cluster_rows)} clusters)')

print(f'\nTop phaG-enriched clusters (baseline {baseline_pct_g:.1f}%):')
print(f'{"cluster_id":34s} {"n":>5s} {"%phaG":>7s} {"top phylum":24s} {"cmh OR":>8s} {"cmh p":>10s}')
for r in cluster_rows[:15]:
    print(f'{r["cluster_id"]:34s} {r["n_genomes"]:>5} {r["pct_phaG"]:>6s}% {r["top_phylum"] + " (" + r["top_phylum_pct"] + "%)":24s} '
          f'{r["cmh_or_phaG"]:>8s} {r["cmh_p_phaG"]:>10s}')

n_sig_elevated = sum(1 for r in cluster_rows if r['cmh_p_phaG'] and float(r['cmh_p_phaG']) < 0.05 and float(r['cmh_or_phaG']) > 1)
n_sig_depleted = sum(1 for r in cluster_rows if r['cmh_p_phaG'] and float(r['cmh_p_phaG']) < 0.05 and float(r['cmh_or_phaG']) < 1)
print(f'\n{n_sig_elevated}/{len(cluster_rows)} clusters significantly ELEVATED for phaG (p<0.05, phylum-controlled)')
print(f'{n_sig_depleted}/{len(cluster_rows)} clusters significantly DEPLETED for phaG (p<0.05, phylum-controlled)')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
STRATEGY_ORDER = ['Neither', 'Beta-oxidation only (J)', 'FAS-linked only (G)', 'Both (G+J)']
STRATEGY_COLORS = {'Neither': '#C9D2CC', 'Beta-oxidation only (J)': '#1E6E7A',
                    'FAS-linked only (G)': '#C2622D', 'Both (G+J)': '#7A9B3E'}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 8.5), dpi=300, gridspec_kw={'width_ratios': [1, 1.1]})

# Panel A: top clusters by size, route composition stacked bar (context, not the headline)
top_by_size = sorted(cluster_rows, key=lambda r: -r['n_genomes'])[:TOP_N_CLUSTERS_PLOT]
top_by_size_sorted = sorted(top_by_size, key=lambda r: r['n_genomes'])
y = np.arange(len(top_by_size_sorted))
bottoms = np.zeros(len(top_by_size_sorted))
counts_map = {
    'Neither': [r['n_neither'] for r in top_by_size_sorted],
    'Beta-oxidation only (J)': [r['n_beta_oxidation_only_J'] for r in top_by_size_sorted],
    'FAS-linked only (G)': [r['n_fas_linked_only_G'] for r in top_by_size_sorted],
    'Both (G+J)': [r['n_both_GJ'] for r in top_by_size_sorted],
}
for strat in STRATEGY_ORDER:
    vals = np.array([100 * v / r['n_genomes'] for v, r in zip(counts_map[strat], top_by_size_sorted)])
    axA.barh(y, vals, left=bottoms, color=STRATEGY_COLORS[strat], height=0.68, zorder=3, label=strat)
    bottoms += vals
axA.set_yticks(y)
axA.set_yticklabels([f'{r["cluster_id"][-12:]}\n{r["top_phylum"]}, n={r["n_genomes"]}' for r in top_by_size_sorted], fontsize=7.8)
axA.set_xlim(0, 100)
axA.set_xlabel('% of single-copy genomes in this cluster')
axA.set_title(f'A. Largest {TOP_N_CLUSTERS_PLOT} clusters (single-copy genomes)', fontsize=12.5, fontweight='bold', loc='left')
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.spines['left'].set_visible(False)
axA.tick_params(axis='y', length=0)
axA.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axA.set_axisbelow(True)

# Panel B: the actual "hunt for pairing" result -- every qualifying cluster's
# %phaG vs. its size, phylum-controlled significance marked
sizes = [r['n_genomes'] for r in cluster_rows]
pct_g = [float(r['pct_phaG']) for r in cluster_rows]
sig_elevated = [r['cmh_p_phaG'] and float(r['cmh_p_phaG']) < 0.05 and float(r['cmh_or_phaG']) > 1 for r in cluster_rows]
sig_depleted = [r['cmh_p_phaG'] and float(r['cmh_p_phaG']) < 0.05 and float(r['cmh_or_phaG']) < 1 for r in cluster_rows]
colors_b = ['#C2622D' if se else ('#3E7CA6' if sd else '#B7C1BC') for se, sd in zip(sig_elevated, sig_depleted)]
axB.scatter(sizes, pct_g, c=colors_b, s=28, alpha=0.85, zorder=3, edgecolor='white', linewidth=0.3)
axB.axhline(baseline_pct_g, color='#5B6E70', linewidth=1.2, linestyle='--', zorder=2)
axB.text(max(sizes) * 0.98, baseline_pct_g, f' baseline: {baseline_pct_g:.1f}%', color='#5B6E70',
          fontsize=9, va='bottom', ha='right', style='italic')
axB.set_xscale('log')
axB.set_xlabel(f'Cluster size (single-copy genomes, log; clusters with >={MIN_CLUSTER_N} shown)')
axB.set_ylabel('% phaG-positive in this cluster')
axB.set_title('B. Every qualifying cluster: %phaG vs. baseline', fontsize=12.5, fontweight='bold', loc='left')
axB.text(0.02, 0.97, f'{n_sig_elevated} clusters significantly elevated (orange)\n{n_sig_depleted} significantly depleted (blue)\np<0.05, phylum-controlled CMH',
          transform=axB.transAxes, ha='left', va='top', fontsize=8.7, color='#3A4A46',
          bbox=dict(boxstyle='round', facecolor='#F2F5F3', edgecolor='#C8D2CD'))
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
axB.grid(True, which='both', color='#E4E8E5', linewidth=0.5, zorder=0)
axB.set_axisbelow(True)

handles = [plt.Rectangle((0, 0), 1, 1, color=STRATEGY_COLORS[s]) for s in STRATEGY_ORDER]
fig.legend(handles, STRATEGY_ORDER, loc='lower center', ncol=4, fontsize=9.5, frameon=False, bbox_to_anchor=(0.28, -0.02))

fig.suptitle('Does phaC cluster identity pair with mcl-PHA precursor route?', fontsize=16, fontweight='bold', y=1.0)
fig.text(0.5, 0.955,
          f'Single-copy genomes only (n={n_final:,}, one phaC cluster per genome, unambiguous pairing) -- '
          f'{n_conflict:,} genomes excluded where the two source pipelines disagreed on copy count.',
          ha='center', fontsize=9.3, color='#5B6E70')
fig.subplots_adjust(left=0.14, right=0.98, top=0.87, bottom=0.14, wspace=0.32)

out_path = OUT / 'phac_cluster_vs_precursor_route.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_cluster_vs_precursor_route.pdf', facecolor='white')
