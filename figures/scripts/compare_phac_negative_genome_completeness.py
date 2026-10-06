"""Is "rich PHA pathway, no phaC" just assembly incompleteness?

The UpSet plot's right-hand group (plot_upset.py) is genomes with more
than 5 PHA families and no synthase found by any method -- including a
whole-proteome mmseqs search. The obvious deflationary explanation is
that those genomes are simply less complete, so phaC was never
assembled. That is testable: every genome in both groups has a CheckM/
anvio completeness in the per-family metadata tables.

Groups compared, defined exactly as the UpSet plot defines them:
  phaC-negative, PHA-rich   no phaC in the uncapped membership list, none
                            among the missed-synthase search's triad-complete
                            hits, and more than 5 other PHA families present
  phaC-positive             a phaC by either of those two routes

Completeness comes from the 15 *_unique_targets_with_metadata.tsv files,
whose `completeness` column is per genome, not per target. One caveat
that cannot be fixed from these files and is reported rather than hidden:
those tables list at most 5 genomes per identical protein, so genomes
hidden by that cap have no completeness row. This costs coverage on the
phaC-positive side specifically (the cap is what hid 7,714 of them), so
the script prints coverage for both groups -- read the comparison only if
both are high.

Usage:
    python figures/scripts/compare_phac_negative_genome_completeness.py
"""
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = ROOT / 'figures/phac_negative_genome_completeness.tsv'

# ------------------------------------------------------------------ groups
bad = _phac_qc.load_bad_targets()
phac_genomes = set()
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        if target_id not in bad:
            phac_genomes.add(genome)
with open(ROOT / 'phac_recovery/novel_phac_candidate_verdicts.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        phac_genomes.add(row['genome'])

pha_rich_no_phac = set()
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        families = {k[2:] for k, v in row.items() if k.startswith('n_pha') and int(v) > 0} - {'phaC'}
        if row['genome'] not in phac_genomes and len(families) > 5:
            pha_rich_no_phac.add(row['genome'])
print(f'{len(phac_genomes):,} phaC-positive genomes; {len(pha_rich_no_phac):,} PHA-rich genomes with no phaC')

# ------------------------------------------------------------- completeness
completeness = {}
for path in sorted(FA.glob('*_unique_targets_with_metadata.tsv')):
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            value = row.get('completeness', '')
            if not value:
                continue
            try:
                completeness[row['genome']] = float(value)
            except ValueError:
                continue
print(f'completeness known for {len(completeness):,} genomes across the 15 family tables')


def summarize(name, genomes):
    values = sorted(completeness[g] for g in genomes if g in completeness)
    coverage = 100 * len(values) / len(genomes)
    q1, q3 = statistics.quantiles(values, n=4)[0], statistics.quantiles(values, n=4)[2]
    print(f'\n{name}')
    print(f'  {len(genomes):>7,} genomes, {len(values):>7,} with completeness ({coverage:.1f}% coverage)')
    print(f'  median {statistics.median(values):.1f}%   IQR {q1:.1f}-{q3:.1f}%   mean {statistics.mean(values):.1f}%')
    for cut in (50, 70, 90):
        print(f'  >={cut}% complete: {100 * sum(v >= cut for v in values) / len(values):.1f}%')
    return values


positive = summarize('phaC-positive', phac_genomes)
negative = summarize('PHA-rich, no phaC', pha_rich_no_phac)


def mann_whitney_u(a, b):
    """Rank-sum with a normal approximation and tie correction -- n is in the
    thousands, so the approximation is sound and this avoids a scipy dependency
    the rest of figures/scripts does not take."""
    combined = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0.0] * len(combined)
    i = 0
    tie_term = 0.0
    while i < len(combined):
        j = i
        while j + 1 < len(combined) and combined[j + 1][0] == combined[i][0]:
            j += 1
        average_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = average_rank
        t = j - i + 1
        tie_term += t ** 3 - t
        i = j + 1
    rank_sum_a = sum(rank for rank, (_, group) in zip(ranks, combined) if group == 0)
    n1, n2 = len(a), len(b)
    u1 = rank_sum_a - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    n = n1 + n2
    sd_u = ((n1 * n2 / 12) * ((n + 1) - tie_term / (n * (n - 1)))) ** 0.5
    z = (u1 - mean_u) / sd_u
    # common-language effect size: P(a random phaC-positive genome is more complete)
    return z, u1 / (n1 * n2)


z, probability = mann_whitney_u(positive, negative)
difference = statistics.median(negative) - statistics.median(positive)
print(f'\nmedian difference (no-phaC minus phaC+): {difference:+.1f} percentage points')
print(f'Mann-Whitney z = {z:.1f}; P(phaC-positive genome is the more complete of a random pair) = {probability:.3f}')

with open(OUT, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['group', 'n_genomes', 'n_with_completeness', 'median', 'q1', 'q3', 'mean',
                'pct_ge_50', 'pct_ge_70', 'pct_ge_90'])
    for name, genomes, values in (('phaC-positive', phac_genomes, positive),
                                  ('pha_rich_no_phaC', pha_rich_no_phac, negative)):
        q1, q3 = statistics.quantiles(values, n=4)[0], statistics.quantiles(values, n=4)[2]
        w.writerow([name, len(genomes), len(values), f'{statistics.median(values):.2f}', f'{q1:.2f}', f'{q3:.2f}',
                    f'{statistics.mean(values):.2f}'] +
                   [f'{100 * sum(v >= cut for v in values) / len(values):.2f}' for cut in (50, 70, 90)])
print(f'\nwrote {OUT}')
