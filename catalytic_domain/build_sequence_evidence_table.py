"""Builds catalytic_domain/phac_sequence_evidence.tsv: one row per verified
phaC candidate target, giving its single best mmseqs2 hit against the
audited reference set, joined to its §3 evidence tier.

No script in the repo previously built this file (it was assembled ad hoc
earlier in the project) -- writing one now because the 2026-09-22 phaC
reference-query fix (PHA_CLEAN_RESULTS.md section 2) makes the old table
doubly stale: not only did the candidate-genome universe shrink (68,424 ->
31,464), the REFERENCE side changed too -- all 152 newly-excluded
accessions were ALSO present in the 1,875-sequence "audited reference
set" itself (confirmed directly: every one of the 152 appears in
structure_prediction/audited_reference_set.faa's headers). That reference
FASTA predates today's fix and is not rebuilt here (deferred along with
section 6/ESMFold, which also consumes it -- see PHA_CLEAN_RESULTS.md).
This script does not need the FASTA rebuilt, though: it works directly
from the raw mmseqs2 hits table (catalytic_domain/audited_search/hits.tsv,
already computed against the full pre-fix reference set) and filters OUT
any hit whose target (reference) accession is one of the 152, which is
equivalent to having searched against a clean reference set in the first
place, just without re-running mmseqs2.

Verified before use, not assumed: of the 33,789 target_ids in the
corrected candidate set, only 9 are entirely absent from hits.tsv's query
universe (101,439 target_ids) -- those 9 genuinely had zero hits against
any of the 1,875 references in the original search below whatever
e-value/sensitivity cutoff was used, not a coverage gap in hits.tsv
itself, so no mmseqs2 re-run is needed.

Best hit per candidate: highest bitscore, ties broken by lowest e-value
(both reported by mmseqs2's own convertalis output, no independent
recomputation).

Usage:
    python catalytic_domain/build_sequence_evidence_table.py
"""
import csv
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'figures/scripts'))
import _phac_qc

# ---------------------------------------------------------------------
# 1. corrected candidate genome/target universe + evidence tier per genome
#    (re-filtered against the CURRENT bad-query list, same pattern as
#    plot_phac_triad_hmm_pathway_groups.py -- never trust a pre-filtered
#    pkl to still be in sync)
# ---------------------------------------------------------------------
genome_targets_raw = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
bad_targets = _phac_qc.load_bad_targets()
genome_targets = {}
target_genome = {}
for g, targets in genome_targets_raw.items():
    kept = targets - bad_targets
    if kept:
        genome_targets[g] = kept
        for t in kept:
            target_genome[t] = g
print(f'{len(genome_targets):,} corrected candidate genomes, {len(target_genome):,} candidate targets')

genome_tier = pickle.load(open('/tmp/phac_verified_triad_hmm_group.pkl', 'rb'))
n_no_tier = sum(1 for g in genome_targets if g not in genome_tier)
if n_no_tier:
    print(f'WARNING: {n_no_tier} candidate genomes have no tier in phac_verified_triad_hmm_group.pkl -- '
          f'rerun plot_phac_triad_hmm_pathway_groups.py first if this is nonzero.')

# ---------------------------------------------------------------------
# 2. bad reference (target-side) accessions -- the 219 currently excluded,
#    same list bad_targets was resolved from, but keyed by ACCESSION here
#    since hits.tsv's target column is a bare UniProt accession, not a
#    "UNIPROT:xxx"-prefixed best_query string.
# ---------------------------------------------------------------------
bad_ref_accs = {q.replace('UNIPROT:', '') for q in _phac_qc.BAD_QUERIES}

# ---------------------------------------------------------------------
# 3. best hit per candidate, from the existing raw mmseqs2 hits table,
#    filtered on both sides (query = corrected candidate set, target =
#    corrected reference set)
# ---------------------------------------------------------------------
best = {}  # target_id -> hit row dict
with open(ROOT / 'catalytic_domain/audited_search/hits.tsv', newline='') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        q = row['query']
        if q not in target_genome:
            continue
        if row['target'] in bad_ref_accs:
            continue
        cur = best.get(q)
        if cur is None or float(row['bits']) > float(cur['bits']) or (
                float(row['bits']) == float(cur['bits']) and float(row['evalue']) < float(cur['evalue'])):
            best[q] = row

print(f'{len(best):,}/{len(target_genome):,} corrected candidates have >=1 surviving hit against the corrected reference set')

# ---------------------------------------------------------------------
# 4. write the table
# ---------------------------------------------------------------------
OUT_PATH = ROOT / 'catalytic_domain/phac_sequence_evidence.tsv'
FIELDNAMES = ['target_id', 'genome', 'evidence_tier', 'protein_length', 'ref_accession', 'ref_length',
              'pident', 'alnlen', 'evalue', 'qcov', 'tcov', 'length_ratio_query_to_ref']

rows_out = []
for target_id, hit in best.items():
    genome = target_genome[target_id]
    tier = genome_tier.get(genome, 'unknown')
    qlen = float(hit['qlen'])
    tlen = float(hit['tlen'])
    rows_out.append({
        'target_id': target_id, 'genome': genome, 'evidence_tier': tier,
        'protein_length': hit['qlen'], 'ref_accession': hit['target'], 'ref_length': hit['tlen'],
        'pident': hit['pident'], 'alnlen': hit['alnlen'], 'evalue': hit['evalue'],
        'qcov': hit['qcov'], 'tcov': hit['tcov'],
        'length_ratio_query_to_ref': f'{qlen / tlen:.4f}' if tlen else '',
    })

with open(OUT_PATH, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)

n_genomes = len(set(r['genome'] for r in rows_out))
print(f'wrote {OUT_PATH} ({len(rows_out):,} rows, {n_genomes:,} distinct genomes)')

from collections import Counter
tier_counts = Counter(r['evidence_tier'] for r in rows_out)
for tier, n in tier_counts.most_common():
    print(f'  {tier}: {n:,}')
