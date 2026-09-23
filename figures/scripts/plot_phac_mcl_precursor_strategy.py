"""mcl-PHA precursor-supply strategy: phaG (de novo fatty-acid-synthesis
route) vs. phaJ (beta-oxidation/scavenging route) vs. the canonical phaA+
phaB (scl-PHA) route, mapped onto full pathway architecture rather than
marginal presence/absence -- and asked whether taxonomy or habitat drives
which strategy a genome carries.

Biology, briefly: phaA+phaB (thiolase + reductase) condenses acetyl-CoA
into (R)-3-hydroxybutyryl-CoA, the canonical scl-PHA (e.g. PHB) precursor
route. For mcl-PHA, a genome needs a DIFFERENT route to generate longer
(R)-3-hydroxyacyl-CoA monomers, and there are two independent ways to get
there: phaJ ((R)-specific enoyl-CoA hydratase) diverts an intermediate out
of beta-oxidation -- i.e. SCAVENGES existing fatty acids, from the
environment or a host -- while phaG (3-hydroxyacyl-ACP:CoA transacylase)
pulls a monomer out of a genome's OWN fatty-acid-synthesis pathway
instead, i.e. builds precursor DE NOVO from acetyl-CoA/malonyl-CoA,
without needing any external fatty acid supply at all. Since scavenging is
normally cheaper than de novo synthesis, phaG being common somewhere would
be a real signal that fatty acids simply are not available there to
scavenge -- a lipid-poor niche.

Restricted to phaC-positive genomes only (68,424 -- this project's
verified scope, PHA_CLEAN_RESULTS.md section 2): phaG/phaJ presence in a
genome with no PHA synthase at all does not describe an actual mcl-PHA
strategy, just incidental fatty-acid-metabolism gene content.

Usage:
    python figures/scripts/plot_phac_mcl_precursor_strategy.py

Outputs:
    figures/phac_mcl_precursor_strategy.png / .pdf
    figures/phac_mcl_strategy_by_phylum.tsv
    figures/phac_mcl_strategy_by_habitat.tsv       (includes phylum-controlled CMH test for phaG)
"""
import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _stats_utils import mantel_haenszel

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

EXCLUDE_HABITATS = {  # keep in sync with plot_phac_pct_by_ocean_habitat.py
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
MIN_HABITAT_N = 200
MIN_PHYLUM_N = 200   # phyla smaller than this are folded into "Other" for panels B/legend readability
TOP_N_PHYLA = 10

# ---------------------------------------------------------------------
# 1. per-genome A/B/G/J presence, phylum -- phaC-positive genomes only
# ---------------------------------------------------------------------
genomes = {}  # genome -> dict(ab, g, j, phylum)
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if int(row['n_phaC']) <= 0:
            continue
        genomes[row['genome']] = {
            'ab': int(row['n_phaA']) > 0 and int(row['n_phaB']) > 0,
            'g': int(row['n_phaG']) > 0,
            'j': int(row['n_phaJ']) > 0,
            'phylum': row['gtdb_phylum'],
        }

n_total = len(genomes)
print(f'{n_total:,} phaC-positive genomes')


def mcl_strategy(rec):
    if rec['g'] and rec['j']:
        return 'Both (G+J)'
    if rec['g']:
        return 'De novo only (G)'
    if rec['j']:
        return 'Scavenging only (J)'
    return 'Neither'


# ---------------------------------------------------------------------
# 2. full A x G x J combination table (panel A)
# ---------------------------------------------------------------------
combo_counts = Counter()
strategy_counts = Counter()
for rec in genomes.values():
    label = ('AB' if rec['ab'] else '') + ('G' if rec['g'] else '') + ('J' if rec['j'] else '')
    combo_counts[label or '(none)'] += 1
    strategy_counts[mcl_strategy(rec)] += 1

print('\nFull A/B/G/J combinations (phaC-positive genomes):')
for label, n in combo_counts.most_common():
    print(f'  {label:10s} n={n:6,d}  ({100*n/n_total:.1f}%)')
print('\nMCL precursor strategy (AB-independent):')
for label, n in sorted(strategy_counts.items(), key=lambda kv: -kv[1]):
    print(f'  {label:22s} n={n:6,d}  ({100*n/n_total:.1f}%)')

# most phaG+ genomes also carry phaJ -- worth surfacing directly, not just implied by the table
n_g = sum(1 for r in genomes.values() if r['g'])
n_g_and_j = sum(1 for r in genomes.values() if r['g'] and r['j'])
print(f'\n{n_g_and_j:,}/{n_g:,} phaG-positive genomes ({100*n_g_and_j/n_g:.1f}%) ALSO carry phaJ -- '
      f'phaG mostly co-occurs with phaJ rather than replacing it.')

# ---------------------------------------------------------------------
# 3. by phylum (panel B) -- genome_family_matrix.tsv has gtdb_phylum for
#    every phaC-positive genome directly, no join/coverage gap here.
# ---------------------------------------------------------------------
phylum_counts = Counter(r['phylum'] for r in genomes.values())
top_phyla = [p for p, _ in phylum_counts.most_common() if phylum_counts[p] >= MIN_PHYLUM_N][:TOP_N_PHYLA]
print(f'\n{len(top_phyla)} phyla with >= {MIN_PHYLUM_N} phaC-positive genomes (of {len(phylum_counts)} total)')

phylum_rows = []
for p in sorted(phylum_counts, key=lambda p: -phylum_counts[p]):
    recs = [r for r in genomes.values() if r['phylum'] == p]
    n = len(recs)
    n_g_p = sum(1 for r in recs if r['g'])
    n_j_p = sum(1 for r in recs if r['j'])
    n_ab_p = sum(1 for r in recs if r['ab'])
    phylum_rows.append({
        'phylum': p, 'n_phac_positive': n,
        'pct_phaG': f'{100*n_g_p/n:.2f}', 'pct_phaJ': f'{100*n_j_p/n:.2f}', 'pct_AB': f'{100*n_ab_p/n:.2f}',
        'n_phaG': n_g_p, 'n_phaJ': n_j_p, 'n_AB': n_ab_p,
    })
phylum_out = OUT / 'phac_mcl_strategy_by_phylum.tsv'
with open(phylum_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(phylum_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(phylum_rows)
print(f'wrote {phylum_out} ({len(phylum_rows)} phyla)')

# ---------------------------------------------------------------------
# 4. by habitat (panel C) + phylum-controlled CMH test for phaG presence
# ---------------------------------------------------------------------
genome_habitat = {}
habitat_totals_all = defaultdict(int)  # among ALL screened genomes (for the MIN_HABITAT_N gate, same convention as habitat figure)
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        c = row['ecosystem_compartment']
        if c in EXCLUDE_HABITATS:
            continue
        habitat_totals_all[c] += 1
        if row['genome'] in genomes:
            genome_habitat[row['genome']] = c

valid_habitats = [c for c, n in habitat_totals_all.items() if n >= MIN_HABITAT_N]
n_with_habitat = sum(1 for g in genomes if g in genome_habitat)
print(f'\n{n_with_habitat:,}/{n_total:,} phaC-positive genomes ({100*n_with_habitat/n_total:.1f}%) have a usable habitat label')

habitat_rows = []
for hab in valid_habitats:
    hab_genomes = [g for g, h in genome_habitat.items() if h == hab]
    if len(hab_genomes) < 20:  # too few phaC-positive genomes in this habitat to report a stable strategy mix
        continue
    recs = [genomes[g] for g in hab_genomes]
    n = len(recs)
    n_g_h = sum(1 for r in recs if r['g'])
    n_j_h = sum(1 for r in recs if r['j'])
    strat = Counter(mcl_strategy(r) for r in recs)

    # phylum-controlled CMH test: is phaG presence still associated with
    # being in this habitat (vs. the rest of the phaC-positive, habitat-
    # labeled dataset) once phylum is controlled for -- same machinery,
    # same reasoning, as phac_habitat_phylum_controlled_test.py's phaC-vs-
    # habitat test, just with the trait swapped from "is phaC positive" to
    # "is phaG positive" (already restricted to phaC-positive genomes).
    phylum_in = defaultdict(lambda: [0, 0])
    phylum_out = defaultdict(lambda: [0, 0])
    for g, h in genome_habitat.items():
        rec = genomes[g]
        bucket = phylum_in if h == hab else phylum_out
        bucket[rec['phylum']][0 if rec['g'] else 1] += 1
    strata = []
    for p in set(phylum_in) | set(phylum_out):
        ia, ib = phylum_in.get(p, [0, 0])
        oa, ob = phylum_out.get(p, [0, 0])
        strata.append((ia, ib, oa, ob))
    or_mh, cmh_stat, cmh_p, n_strata = mantel_haenszel(strata)

    a = n_g_h
    b = n - n_g_h
    c_ = n_g - n_g_h
    d = (n_with_habitat - n) - c_
    raw_or, raw_p = fisher_exact([[a, b], [c_, d]], alternative='two-sided')

    habitat_rows.append({
        'habitat': hab, 'n_phac_positive': n,
        'pct_phaG': f'{100*n_g_h/n:.2f}', 'pct_phaJ': f'{100*n_j_h/n:.2f}',
        'n_neither': strat.get('Neither', 0), 'n_scavenging_only_J': strat.get('Scavenging only (J)', 0),
        'n_de_novo_only_G': strat.get('De novo only (G)', 0), 'n_both_GJ': strat.get('Both (G+J)', 0),
        'raw_fisher_or_phaG': f'{raw_or:.3f}', 'raw_fisher_p_phaG': f'{raw_p:.3g}',
        'cmh_or_phaG': f'{or_mh:.3f}' if or_mh is not None else '',
        'cmh_p_phaG': f'{cmh_p:.3g}' if cmh_p is not None else '',
        'n_strata_used': n_strata,
    })

habitat_rows.sort(key=lambda r: -float(r['pct_phaG']))
habitat_out = OUT / 'phac_mcl_strategy_by_habitat.tsv'
with open(habitat_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(habitat_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(habitat_rows)
print(f'wrote {habitat_out} ({len(habitat_rows)} habitats)')
print(f'\n{"habitat":32s} {"n":>7s} {"%phaG":>7s} {"%phaJ":>7s} {"cmh OR":>8s} {"cmh p":>10s}')
for r in habitat_rows:
    print(f'{r["habitat"]:32s} {r["n_phac_positive"]:>7} {r["pct_phaG"]:>6s}% {r["pct_phaJ"]:>6s}% '
          f'{r["cmh_or_phaG"]:>8s} {r["cmh_p_phaG"]:>10s}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
STRATEGY_ORDER = ['Neither', 'Scavenging only (J)', 'De novo only (G)', 'Both (G+J)']
STRATEGY_COLORS = {'Neither': '#C9D2CC', 'Scavenging only (J)': '#1E6E7A',
                    'De novo only (G)': '#C2622D', 'Both (G+J)': '#7A9B3E'}

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig = plt.figure(figsize=(15, 13.5), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.7], hspace=0.3, wspace=0.32)
axA = fig.add_subplot(gs[0, :])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])

# Panel A: full combo bar chart, log y
combo_order = sorted(combo_counts, key=lambda k: -combo_counts[k])
axA.bar(combo_order, [combo_counts[k] for k in combo_order], color='#3A6B63', zorder=3)
axA.set_yscale('log')
for i, k in enumerate(combo_order):
    axA.text(i, combo_counts[k] * 1.15, f'{combo_counts[k]:,}', ha='center', fontsize=8.7, color='#20302C')
axA.set_ylabel('Genomes (log)')
axA.set_title('A. Full pathway architecture: canonical (AB) x scavenging (J) x de novo (G)',
               fontsize=12.5, fontweight='bold', loc='left')
axA.text(0.99, 0.95, f'n={n_total:,} phaC-positive genomes\n{n_g_and_j:,}/{n_g:,} phaG+ genomes ({100*n_g_and_j/n_g:.0f}%) also carry phaJ',
          transform=axA.transAxes, ha='right', va='top', fontsize=9, color='#3A4A46',
          bbox=dict(boxstyle='round', facecolor='#F2F5F3', edgecolor='#C8D2CD'))
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.grid(axis='y', which='both', color='#E4E8E5', linewidth=0.5, zorder=0)
axA.set_axisbelow(True)


def stacked_strategy_bar(ax, categories, counts_by_cat, title, label_fontsize=9.5):
    n_cat = len(categories)
    y = np.arange(n_cat)  # explicit numeric positions -- repeated barh() calls against the
    # same string-category list (once per stacked segment) misaligned tick labels vs bars
    # in testing; numeric positions + explicit set_yticklabels sidesteps that entirely.
    bottoms = np.zeros(n_cat)
    for strat in STRATEGY_ORDER:
        vals = np.array([100 * counts_by_cat[cat].get(strat, 0) / sum(counts_by_cat[cat].values()) for cat in categories])
        ax.barh(y, vals, left=bottoms, color=STRATEGY_COLORS[strat], height=0.68, zorder=3, label=strat)
        bottoms += vals
    ax.set_yticks(y)
    ax.set_yticklabels(categories)
    ax.set_ylim(-0.6, n_cat - 0.4)
    ax.set_xlim(0, 100)
    ax.set_xlabel('% of phaC-positive genomes')
    ax.set_title(title, fontsize=12.5, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, labelsize=label_fontsize)
    ax.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


# Panel B: by phylum
phylum_counts_by_cat = {}
for p in top_phyla:
    recs = [r for r in genomes.values() if r['phylum'] == p]
    phylum_counts_by_cat[p] = Counter(mcl_strategy(r) for r in recs)
phyla_sorted = sorted(top_phyla, key=lambda p: sum(phylum_counts_by_cat[p].values()))
labels_b = [f'{p} (n={phylum_counts[p]:,})' for p in phyla_sorted]
stacked_strategy_bar(axB, labels_b, {labels_b[i]: phylum_counts_by_cat[p] for i, p in enumerate(phyla_sorted)},
                      f'B. MCL precursor strategy by phylum (top {len(top_phyla)})')

# Panel C: by habitat
habitat_counts_by_cat = {}
for r in habitat_rows:
    hab = r['habitat']
    habitat_counts_by_cat[hab] = {
        'Neither': r['n_neither'], 'Scavenging only (J)': r['n_scavenging_only_J'],
        'De novo only (G)': r['n_de_novo_only_G'], 'Both (G+J)': r['n_both_GJ'],
    }
habitats_sorted = sorted(habitat_counts_by_cat, key=lambda h: sum(habitat_counts_by_cat[h].values()))
labels_c = [f'{h} (n={sum(habitat_counts_by_cat[h].values()):,})' for h in habitats_sorted]
stacked_strategy_bar(axC, labels_c, {labels_c[i]: habitat_counts_by_cat[h] for i, h in enumerate(habitats_sorted)},
                      'C. MCL precursor strategy by habitat', label_fontsize=8.5)

handles = [plt.Rectangle((0, 0), 1, 1, color=STRATEGY_COLORS[s]) for s in STRATEGY_ORDER]
fig.legend(handles, STRATEGY_ORDER, loc='lower center', ncol=4, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, -0.01))

fig.suptitle('mcl-PHA precursor-supply strategy: scavenging (phaJ) vs. de novo synthesis (phaG)',
             fontsize=16, fontweight='bold', y=0.995)
fig.text(0.5, 0.965,
          f'phaC-positive genomes only (n={n_total:,}). "De novo" and "Both" bars are phaG-positive; phaG is rare '
          f'overall ({100*n_g/n_total:.1f}% of phaC-positive genomes) and mostly co-occurs with phaJ rather than replacing it.',
          ha='center', fontsize=9.3, color='#5B6E70')
fig.subplots_adjust(left=0.16, right=0.98, top=0.94, bottom=0.06)

out_path = OUT / 'phac_mcl_precursor_strategy.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_mcl_precursor_strategy.pdf', facecolor='white')
