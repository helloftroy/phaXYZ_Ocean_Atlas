"""How many phaC copies does a genome carry, and is that (taxonomically or
habitat-wise) driven by anything real, or mostly assembly/binning noise?

Uses `n_phaC` from genome_family_matrix.tsv directly -- the same column
already authoritative for this project's whole "verified family scope"
(PHA_CLEAN_RESULTS.md section 2), rather than reconstructing a copy count
from the raw NR100-clustered target_id->genome join
(phaC_all_genomes_from_nr100_clusters.tsv): that file's target_id is a
100%-identity-clustered REPRESENTATIVE id shared across every genome
carrying an identical sequence, so counting its rows per genome conflates
"this genome has multiple different phaC genes" with "this exact phaC
sequence happens to also appear in the NR100 catalogue attached to other
genomes' redundant entries" -- not a real second signal. genome_family_matrix.tsv's
n_phaC is a straightforward per-genome gene count and doesn't have this
confound.

Multi-copy phaC is genuinely interesting biology when real (e.g. plasmid-
or duplication-borne extra copies, sometimes with distinct substrate
specificities) but a genome assembled from a metagenome with very high
apparent copy number (this dataset's max is 30) is at least as likely to
be a binning/contamination artifact -- two or more organisms' contigs
merged into one "genome" bin. Checked directly below via each genome's own
CheckM-style contamination score (from phaC_unique_targets_with_metadata_depth.tsv)
rather than assumed.

Rebuilt again 2026-09-24 with two more corrections layered on top of
n_phaC, following section 9.12's structural resolution of the
no_hmm_triad_support population (0/19 have a real catalytic triad by any
mechanism checked -- not phaC): (1) decrement n_phaC by 1 for every
genome carrying one of those 19 directly-excluded target_ids (via
phaC_all_genomes_from_nr100_clusters.tsv's genome join -- 23 genomes
affected, since a few of the 19 are NR100 100%-identity representatives
shared across more than one genome); (2) drop the 85 genomes
figures/phac_multicopy_legitimacy_audit.tsv's composite audit (section
9.3) called `likely_artifact` outright, rather than trusting their
n_phaC at all -- that call already means the genome's own copies mostly
fail the triad-completeness/distinct-cluster checks, i.e. look like
assembly fragmentation, not real paralogs. Both are small corrections
(23 + 85 genomes out of ~31,000) but the right thing to apply now that
section 9.12 gives a concrete, checked reason for the first one instead
of a suspicion.

Usage:
    python figures/scripts/plot_phac_multicopy_genomes.py

Outputs:
    figures/phac_multicopy_genomes.png / .pdf
    figures/phac_multicopy_by_phylum.tsv
    figures/phac_multicopy_by_habitat.tsv    (includes phylum-controlled CMH test for multi-copy)
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
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / 'figures'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

EXCLUDE_HABITATS = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}
MIN_HABITAT_N = 200
MIN_PHYLUM_N = 200
TOP_N_PHYLA = 12

# ---------------------------------------------------------------------
# 1. n_phaC per genome (authoritative) + phylum
# ---------------------------------------------------------------------
genomes = {}  # genome -> dict(n_phac, phylum)
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genomes[row['genome']] = {'n_phac': n, 'phylum': row['gtdb_phylum']}

n_before_corrections = len(genomes)

# correction 1: decrement n_phaC for genomes carrying one of the 19
# section-9.12-excluded target_ids (structurally confirmed to have no
# real catalytic triad, canonical or alternative -- not phaC)
bad_target_hits = defaultdict(int)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in _phac_qc.BAD_TARGET_IDS and row['genome'] in genomes:
            bad_target_hits[row['genome']] += 1
for g, n_bad in bad_target_hits.items():
    genomes[g]['n_phac'] = max(0, genomes[g]['n_phac'] - n_bad)
n_dropped_to_zero = sum(1 for g in bad_target_hits if genomes[g]['n_phac'] == 0)
genomes = {g: r for g, r in genomes.items() if r['n_phac'] > 0}
print(f'correction 1: {len(bad_target_hits)} genomes had >=1 of the 19 excluded target_ids '
      f'(n_phaC decremented accordingly; {n_dropped_to_zero} dropped out of the phaC-positive set entirely)')

# correction 2: drop genomes the section-9.3 composite legitimacy audit
# called likely_artifact outright (assembly-fragmentation signature, not
# trusted at all rather than just decremented)
likely_artifact = set()
with open(OUT / 'phac_multicopy_legitimacy_audit.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['legitimacy_call'] == 'likely_artifact':
            likely_artifact.add(row['genome'])
n_artifact_removed = sum(1 for g in genomes if g in likely_artifact)
genomes = {g: r for g, r in genomes.items() if g not in likely_artifact}
print(f'correction 2: {n_artifact_removed} likely_artifact genomes dropped entirely '
      f'(of {len(likely_artifact)} total in the legitimacy audit)')
print(f'{n_before_corrections:,} -> {len(genomes):,} phaC-positive genomes after both corrections\n')

n_total = len(genomes)
copy_dist = Counter(r['n_phac'] for r in genomes.values())
n_multi = sum(1 for r in genomes.values() if r['n_phac'] >= 2)
print(f'{n_total:,} phaC-positive genomes; {n_multi:,} ({100*n_multi/n_total:.1f}%) carry >=2 phaC copies')
for k in sorted(copy_dist)[:12]:
    print(f'  n_phaC={k:2d}: {copy_dist[k]:6,d}  ({100*copy_dist[k]/n_total:.2f}%)')
n_high = sum(v for k, v in copy_dist.items() if k >= 10)
print(f'  n_phaC>=10: {n_high:,} ({100*n_high/n_total:.2f}%) -- treat as likely binning/contamination-affected, not real 10+ paralogs')

# ---------------------------------------------------------------------
# 2. contamination sanity check (CheckM-style % from phaC_unique_targets_with_metadata_depth.tsv,
#    genome-level -- take first non-empty value seen per genome, contamination
#    is a property of the genome/bin, not of the specific target row)
# ---------------------------------------------------------------------
genome_contam = {}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row['genome']
        if g not in genome_contam and row['contamination']:
            try:
                genome_contam[g] = float(row['contamination'])
            except ValueError:
                pass

BUCKETS = [(1, 1, '1'), (2, 2, '2'), (3, 4, '3-4'), (5, 9, '5-9'), (10, 999, '10+')]
contam_by_bucket = defaultdict(list)
n_with_contam = 0
for g, rec in genomes.items():
    c = genome_contam.get(g)
    if c is None:
        continue
    n_with_contam += 1
    for lo, hi, label in BUCKETS:
        if lo <= rec['n_phac'] <= hi:
            contam_by_bucket[label].append(c)
            break

print(f'\n{n_with_contam:,}/{n_total:,} genomes ({100*n_with_contam/n_total:.1f}%) have a contamination score; mean contamination by copy-count bucket:')
for _, _, label in BUCKETS:
    vals = contam_by_bucket[label]
    if vals:
        print(f'  {label:5s} copies: mean contamination = {sum(vals)/len(vals):.2f}%  (n={len(vals):,})')

# ---------------------------------------------------------------------
# 3. by phylum
# ---------------------------------------------------------------------
phylum_counts = Counter(r['phylum'] for r in genomes.values())
top_phyla = [p for p, _ in phylum_counts.most_common() if phylum_counts[p] >= MIN_PHYLUM_N][:TOP_N_PHYLA]

phylum_rows = []
for p in sorted(phylum_counts, key=lambda p: -phylum_counts[p]):
    recs = [r for r in genomes.values() if r['phylum'] == p]
    n = len(recs)
    n_multi_p = sum(1 for r in recs if r['n_phac'] >= 2)
    mean_copies = sum(r['n_phac'] for r in recs) / n
    phylum_rows.append({'phylum': p, 'n_phac_positive': n, 'pct_multicopy': f'{100*n_multi_p/n:.2f}',
                         'mean_n_phac': f'{mean_copies:.2f}', 'n_multicopy': n_multi_p})
phylum_out = OUT / 'phac_multicopy_by_phylum.tsv'
with open(phylum_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(phylum_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(phylum_rows)
print(f'\nwrote {phylum_out} ({len(phylum_rows)} phyla)')

# ---------------------------------------------------------------------
# 4. by habitat + phylum-controlled CMH test for "multi-copy" (n_phaC>=2)
# ---------------------------------------------------------------------
genome_habitat = {}
habitat_totals_all = defaultdict(int)
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
overall_pct_multi = 100 * n_multi / n_total

habitat_rows = []
for hab in valid_habitats:
    hab_genomes = [g for g, h in genome_habitat.items() if h == hab]
    if len(hab_genomes) < 20:
        continue
    recs = [genomes[g] for g in hab_genomes]
    n = len(recs)
    n_multi_h = sum(1 for r in recs if r['n_phac'] >= 2)

    phylum_in = defaultdict(lambda: [0, 0])
    phylum_out = defaultdict(lambda: [0, 0])
    for g, h in genome_habitat.items():
        rec = genomes[g]
        bucket = phylum_in if h == hab else phylum_out
        bucket[rec['phylum']][0 if rec['n_phac'] >= 2 else 1] += 1
    strata = []
    for p in set(phylum_in) | set(phylum_out):
        ia, ib = phylum_in.get(p, [0, 0])
        oa, ob = phylum_out.get(p, [0, 0])
        strata.append((ia, ib, oa, ob))
    or_mh, cmh_stat, cmh_p, n_strata = mantel_haenszel(strata)

    a = n_multi_h
    b = n - n_multi_h
    c_ = n_multi - n_multi_h
    d = (n_with_habitat - n) - c_
    raw_or, raw_p = fisher_exact([[a, b], [c_, d]], alternative='two-sided')

    habitat_rows.append({
        'habitat': hab, 'n_phac_positive': n, 'pct_multicopy': f'{100*n_multi_h/n:.2f}',
        'raw_fisher_or': f'{raw_or:.3f}', 'raw_fisher_p': f'{raw_p:.3g}',
        'cmh_or': f'{or_mh:.3f}' if or_mh is not None else '',
        'cmh_p': f'{cmh_p:.3g}' if cmh_p is not None else '',
        'n_strata_used': n_strata,
    })

habitat_rows.sort(key=lambda r: -float(r['pct_multicopy']))
habitat_out = OUT / 'phac_multicopy_by_habitat.tsv'
with open(habitat_out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(habitat_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(habitat_rows)
print(f'wrote {habitat_out} ({len(habitat_rows)} habitats)')
print(f'\noverall %multicopy = {overall_pct_multi:.1f}%')
print(f'{"habitat":32s} {"n":>7s} {"%multi":>7s} {"cmh OR":>8s} {"cmh p":>10s}')
for r in habitat_rows:
    print(f'{r["habitat"]:32s} {r["n_phac_positive"]:>7} {r["pct_multicopy"]:>6s}% {r["cmh_or"]:>8s} {r["cmh_p"]:>10s}')

# ---------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, axes = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
axA, axB, axC, axD = axes.flat

# Panel A: copy-count distribution
copy_labels = sorted(copy_dist)
axA.bar([str(k) for k in copy_labels], [copy_dist[k] for k in copy_labels], color='#3A6B63', zorder=3)
axA.set_yscale('log')
axA.set_xlabel('phaC copies in genome')
axA.set_ylabel('Genomes (log)')
axA.set_title('A. phaC copy-count distribution', fontsize=12.5, fontweight='bold', loc='left')
axA.text(0.97, 0.95, f'n={n_total:,} genomes\n{n_multi:,} ({100*n_multi/n_total:.1f}%) have >=2 copies',
          transform=axA.transAxes, ha='right', va='top', fontsize=9, color='#3A4A46',
          bbox=dict(boxstyle='round', facecolor='#F2F5F3', edgecolor='#C8D2CD'))
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.grid(axis='y', which='both', color='#E4E8E5', linewidth=0.5, zorder=0)
axA.set_axisbelow(True)
axA.set_xticks(range(0, len(copy_labels), 2))
axA.set_xticklabels([str(copy_labels[i]) for i in range(0, len(copy_labels), 2)])

# Panel B: mean contamination by copy-count bucket
bucket_labels = [b[2] for b in BUCKETS]
means = [sum(contam_by_bucket[l]) / len(contam_by_bucket[l]) if contam_by_bucket[l] else 0 for l in bucket_labels]
ns = [len(contam_by_bucket[l]) for l in bucket_labels]
sems = [np.std(contam_by_bucket[l], ddof=1) / np.sqrt(len(contam_by_bucket[l])) if len(contam_by_bucket[l]) > 1 else 0 for l in bucket_labels]
axB.bar(bucket_labels, means, yerr=sems, color='#9E3B3B', alpha=0.85, capsize=4, zorder=3)
axB.set_ylim(0, max(means) * 1.8)
for i, (m, n) in enumerate(zip(means, ns)):
    axB.text(i, m + 0.22, f'n={n:,}', ha='center', fontsize=8.5, color='#3A4A46')
axB.set_xlabel('phaC copies in genome')
axB.set_ylabel('Mean contamination (%, CheckM-style)')
axB.set_title('B. Does copy count track with contamination?', fontsize=12.5, fontweight='bold', loc='left')
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
axB.grid(axis='y', color='#E4E8E5', linewidth=0.6, zorder=0)
axB.set_axisbelow(True)

# Panel C: %multicopy by phylum
phyla_for_plot = sorted(top_phyla, key=lambda p: sum(1 for r in genomes.values() if r['phylum'] == p and r['n_phac'] >= 2) /
                          max(1, sum(1 for r in genomes.values() if r['phylum'] == p)))
pct_by_phylum = []
for p in phyla_for_plot:
    recs = [r for r in genomes.values() if r['phylum'] == p]
    pct_by_phylum.append(100 * sum(1 for r in recs if r['n_phac'] >= 2) / len(recs))
colors_c = ['#1E6E7A' if v >= overall_pct_multi else '#C9622D' for v in pct_by_phylum]
axC.barh([f'{p} (n={phylum_counts[p]:,})' for p in phyla_for_plot], pct_by_phylum, color=colors_c, height=0.68, zorder=3)
axC.axvline(overall_pct_multi, color='#5B6E70', linewidth=1.2, linestyle='--', zorder=2)
axC.set_xlabel('% of phaC-positive genomes with >=2 copies')
axC.set_title(f'C. Multi-copy prevalence by phylum (top {len(top_phyla)})', fontsize=12.5, fontweight='bold', loc='left')
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)
axC.spines['left'].set_visible(False)
axC.tick_params(axis='y', length=0, labelsize=9.5)
axC.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axC.set_axisbelow(True)

# Panel D: %multicopy by habitat
habs_for_plot = sorted(habitat_rows, key=lambda r: float(r['pct_multicopy']))
labels_d = [f'{r["habitat"]} (n={r["n_phac_positive"]:,})' for r in habs_for_plot]
vals_d = [float(r['pct_multicopy']) for r in habs_for_plot]
colors_d = ['#1E6E7A' if v >= overall_pct_multi else '#C9622D' for v in vals_d]
axD.barh(labels_d, vals_d, color=colors_d, height=0.68, zorder=3)
axD.axvline(overall_pct_multi, color='#5B6E70', linewidth=1.2, linestyle='--', zorder=2)
axD.text(overall_pct_multi, len(labels_d) - 0.3, f'  overall: {overall_pct_multi:.1f}%', color='#5B6E70',
          fontsize=8.7, va='center', ha='left', style='italic')
axD.set_xlabel('% of phaC-positive genomes with >=2 copies')
axD.set_title('D. Multi-copy prevalence by habitat', fontsize=12.5, fontweight='bold', loc='left')
axD.spines['top'].set_visible(False)
axD.spines['right'].set_visible(False)
axD.spines['left'].set_visible(False)
axD.tick_params(axis='y', length=0, labelsize=8.5)
axD.grid(axis='x', color='#E4E8E5', linewidth=0.6, zorder=0)
axD.set_axisbelow(True)

for ax in axes.flat:
    pass

fig.suptitle('Multiple phaC copies per genome: how common, and driven by what', fontsize=16, fontweight='bold', y=0.995)
fig.text(0.5, 0.965,
          f'{n_total:,} phaC-positive genomes ({100*n_multi/n_total:.1f}% carry >=2 copies). Panel D bars colored teal/orange by '
          f'above/below the {overall_pct_multi:.1f}% overall rate; see figures/phac_multicopy_by_habitat.tsv for the phylum-controlled CMH test.',
          ha='center', fontsize=9.3, color='#5B6E70')
fig.subplots_adjust(left=0.18, right=0.98, top=0.92, bottom=0.06, hspace=0.32, wspace=0.35)

out_path = OUT / 'phac_multicopy_genomes.png'
fig.savefig(out_path, dpi=300, facecolor='white')
print('saved', out_path)
fig.savefig(OUT / 'phac_multicopy_genomes.pdf', facecolor='white')
