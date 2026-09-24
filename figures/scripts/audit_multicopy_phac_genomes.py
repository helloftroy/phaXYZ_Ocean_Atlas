"""Scriptable version of the manual legitimacy checks used in section 9.1's
high-copy deep dive and section 9.2's Legionella deep dive, run against
EVERY genome with >2 phaC copies (not just the >=5-copy tail or one
genus) -- so the approximate fraction of "real paralogs" vs. "assembly/
fragmentation artifact" genomes can be judged quickly across the whole
multi-copy population, before picking specific cases for a deep-dive
figure.

Per genome, five signals, each independently computed from data already
in this repo (no network calls -- see section 9.2 for why an NCBI check
does not scale past a handful of hand-picked genomes):

  1. pct_triad_complete -- of this genome's own phaC target_ids (via the
     corrected NR100 target_id<->genome join, filtered against the
     current 219-accession bad-target list -- never trust n_phaC from
     genome_family_matrix.tsv alone for this: section 9.2 found genomes
     where it undercounts against this join by 5-6x), what fraction have
     triad_complete=True in catalytic_domain/phac_catalytic_triad.tsv's
     own per-target column. This is the same real per-target check
     section 9.1 originally got wrong by reading it off a genome-level-
     broadcast column instead -- see that section's correction note.

  2. pct_full_length -- fraction of own targets with protein_length_aa
     >= FULL_LENGTH_AA (500aa). Chosen against this project's own
     reference-set statistics (section 4): median catalytic-triad-complete
     candidate length is 563aa, median reference length 591aa -- 500aa is
     a deliberately generous floor below the typical full protein, meant
     to catch obvious fragments (the Legionella W1046 fragments were
     96-367aa) without penalizing genuine natural length variation.

  3. pct_distinct_clusters -- n_distinct_70%-clusters / n_targets among a
     genome's own copies. 1.0 = every copy is a genuinely distinct
     paralog; well below 1.0 = some copies are near-duplicates of each
     other, the pattern section 9.2 found in the fragmented W1046 genome
     (5 distinct clusters for 9 copies) vs. the clean FDAARGOS_200 genome
     (9 distinct clusters for 9 copies).

  4. checkm_completeness / checkm_contamination -- from
     phaC_unique_targets_with_metadata_depth.tsv, available locally for
     both MAG and isolate genomes. Caveat, confirmed directly: this
     project's own completeness figure does not always agree with NCBI's
     own per-assembly CheckM figure -- W1046 (RefSeq-suppressed, NCBI
     completeness 66.75%) shows 96.0% here. Treat this signal as weaker
     than 1-3, not a substitute for them.

  5. n_phaC_contigs / max_phaC_on_single_contig -- among the subset of a
     genome's own targets whose phaC_cluster_sequences.faa "rep=" example
     happens to be this same genome (a target shared identically with
     another genome only records one arbitrary rep, so this is a LOWER
     bound on true contig count/max-per-contig, not exhaustive --
     n_phaC_contigs_known_of reports how many of the genome's targets this
     could even be checked for, so the reader can judge how much to trust
     it per genome): n_phaC_contigs is the number of distinct scaffolds
     those targets fall on; max_phaC_on_single_contig is the largest
     number of them found on any one scaffold. Both matter for reading
     section 9.4's contig-map style finding at scale: a genome with many
     copies spread across many different contigs (high n_phaC_contigs,
     low max_phaC_on_single_contig) is a different biological story than
     one with several copies piled on the same contig (which could be a
     real tandem cluster, like section 9.4's phaE-phaC-phaJ operon, or a
     single mis-assembled/fragmented region -- this table cannot tell
     those apart by itself, only flag genomes worth a closer look).
     This pass is local-only, same as 1-4 (no network calls, run against
     all 3,176 genomes) -- for a genome with incomplete
     n_phaC_contigs_known_of coverage picked for an actual deep dive, a
     direct fetch of that genome's own gene calls (as section 9.4 did for
     CARD22-1) resolves the gap completely; that is not done here for all
     3,176 genomes, only worth doing by hand for whichever few are
     shortlisted next.

Composite call (LEGIT_TRIAD_MIN / LEGIT_CLUSTER_MIN / ARTIFACT_TRIAD_MAX /
ARTIFACT_CLUSTER_MAX below): "likely_legit" if pct_triad_complete >= 0.6
AND pct_distinct_clusters >= 0.6 AND contamination (if known) <= 10%;
"likely_artifact" if pct_triad_complete < 0.3 OR pct_distinct_clusters <
0.3 OR contamination (if known) > 10%; else "uncertain". These thresholds
are a heuristic, deliberately simple and stated here rather than a
black-box score -- calibrated against section 9.2's two hand-checked
extremes (W1046: 0.22 triad / 0.56 clusters -> uncertain-leaning-artifact;
FDAARGOS_200: 0.78 triad / 1.0 clusters -> likely_legit), not fit to the
full dataset. Read the individual columns, not just the call, before
trusting any one genome's classification.

Usage:
    python figures/scripts/audit_multicopy_phac_genomes.py

Outputs:
    figures/phac_multicopy_legitimacy_audit.tsv   (one row per genome, n_phaC>=3)
No figures -- this is a scanning/triage step; specific genomes worth a
deep-dive figure get picked from this table's survivors by hand.
"""
import csv
import sys
from collections import defaultdict, Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _phac_qc

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'

MIN_COPIES = 3          # ">2 phaC", per the request
FULL_LENGTH_AA = 500
LEGIT_TRIAD_MIN = 0.6
LEGIT_CLUSTER_MIN = 0.6
ARTIFACT_TRIAD_MAX = 0.3
ARTIFACT_CLUSTER_MAX = 0.3
CONTAM_FAIL_PCT = 10.0

# ---------------------------------------------------------------------
# 1. candidate genomes: n_phaC >= MIN_COPIES per genome_family_matrix.tsv
# ---------------------------------------------------------------------
genomes = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        n = int(row['n_phaC'])
        if n <= 0:
            continue
        genomes[row['genome']] = row

multi = {g: r for g, r in genomes.items() if int(r['n_phaC']) >= MIN_COPIES}
multi_set = set(multi)
print(f'{len(genomes):,} phaC-positive genomes; {len(multi):,} ({100*len(multi)/len(genomes):.2f}%) have n_phaC>={MIN_COPIES}')

# ---------------------------------------------------------------------
# 2. corrected per-target join (the authoritative copy count -- see docstring)
# ---------------------------------------------------------------------
bad_targets = _phac_qc.load_bad_targets()
genome_targets = defaultdict(set)
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        t, g = row['target_id'], row['genome']
        if g in multi_set and t not in bad_targets:
            genome_targets[g].add(t)

n_no_targets = sum(1 for g in multi if not genome_targets.get(g))
if n_no_targets:
    print(f'WARNING: {n_no_targets} n_phaC>={MIN_COPIES} genomes have zero surviving targets in the corrected join -- excluded below')

all_targets = set()
for s in genome_targets.values():
    all_targets |= s
print(f'{len(all_targets):,} distinct target_ids across these genomes')

# ---------------------------------------------------------------------
# 3. per-target: triad_complete, protein length, cluster0.7 assignment, scaffold
# ---------------------------------------------------------------------
triad_complete = {}
with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in all_targets:
            triad_complete[row['target_id']] = row['triad_complete'] == 'True'

length_aa = {}
with open(ROOT / 'catalytic_domain/phac_sequence_evidence.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in all_targets:
            length_aa[row['target_id']] = int(row['protein_length'])

member_to_cluster = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if row[1] in all_targets:
            member_to_cluster[row[1]] = row[0]

scaffold_of = {}  # target_id -> (genome, scaffold) only when rep IS that genome
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if not line.startswith('>'):
            continue
        parts = line[1:].strip().split(' ', 1)
        tid = parts[0]
        if tid not in all_targets or len(parts) < 2:
            continue
        rep_body = parts[1].split('rep=')[-1].split(';')[0]
        # rep_body = "GENOME-scaffold_N_M" -- genome names themselves can
        # contain hyphens, so split on the LAST "-scaffold_" marker
        if '-scaffold_' not in rep_body:
            continue
        rep_genome, tail = rep_body.rsplit('-scaffold_', 1)
        scaffold_of[tid] = (rep_genome, 'scaffold_' + tail.rsplit('_', 1)[0])

print(f'triad table coverage: {len(triad_complete):,}/{len(all_targets):,}; '
      f'length coverage: {len(length_aa):,}/{len(all_targets):,}; '
      f'cluster coverage: {len(member_to_cluster):,}/{len(all_targets):,}; '
      f'scaffold (own-genome-only) coverage: {sum(1 for t,(g,s) in scaffold_of.items())}/{len(all_targets):,}')

# ---------------------------------------------------------------------
# 4. genome-level quality (completeness/contamination), MAG or isolate alike
# ---------------------------------------------------------------------
genome_qc = {}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        g = row['genome']
        if g in multi_set and g not in genome_qc and row['completeness']:
            genome_qc[g] = (row['completeness'], row['contamination'])

# ---------------------------------------------------------------------
# 5. assemble per-genome rows + composite call
# ---------------------------------------------------------------------
rows_out = []
for g, r in multi.items():
    tids = genome_targets.get(g)
    if not tids:
        continue
    n = len(tids)
    n_triad_known = sum(1 for t in tids if t in triad_complete)
    n_triad_true = sum(1 for t in tids if triad_complete.get(t))
    pct_triad = n_triad_true / n_triad_known if n_triad_known else None

    lens = [length_aa[t] for t in tids if t in length_aa]
    pct_full = sum(1 for l in lens if l >= FULL_LENGTH_AA) / len(lens) if lens else None
    median_len = sorted(lens)[len(lens) // 2] if lens else None

    clusters = set(member_to_cluster.get(t) for t in tids if t in member_to_cluster)
    n_clusters_known = sum(1 for t in tids if t in member_to_cluster)
    pct_distinct = len(clusters) / n_clusters_known if n_clusters_known else None

    own_scaffold_counts = Counter()
    n_scaffold_known = 0
    for t in tids:
        rg_s = scaffold_of.get(t)
        if rg_s and rg_s[0] == g:
            own_scaffold_counts[rg_s[1]] += 1
            n_scaffold_known += 1
    n_phac_contigs = len(own_scaffold_counts)
    max_on_one_contig = max(own_scaffold_counts.values()) if own_scaffold_counts else 0

    comp, contam = genome_qc.get(g, ('', ''))
    contam_f = float(contam) if contam else None

    # composite call -- see module docstring for the exact rule and its calibration
    if pct_triad is None or pct_distinct is None:
        call = 'insufficient_data'
    elif (pct_triad >= LEGIT_TRIAD_MIN and pct_distinct >= LEGIT_CLUSTER_MIN
          and (contam_f is None or contam_f <= CONTAM_FAIL_PCT)):
        call = 'likely_legit'
    elif (pct_triad < ARTIFACT_TRIAD_MAX or pct_distinct < ARTIFACT_CLUSTER_MAX
          or (contam_f is not None and contam_f > CONTAM_FAIL_PCT)):
        call = 'likely_artifact'
    else:
        call = 'uncertain'

    rows_out.append({
        'genome': g, 'gtdb_species': r['gtdb_species'], 'gtdb_genus': r['gtdb_genus'],
        'study_id': r['study_id'], 'ecosystem_compartment': r['ecosystem_compartment'],
        'is_mag': 'MAG' in g,
        'n_phaC_matrix': r['n_phaC'], 'n_targets_corrected': n,
        'pct_triad_complete': f'{pct_triad:.3f}' if pct_triad is not None else '',
        'pct_full_length': f'{pct_full:.3f}' if pct_full is not None else '',
        'median_protein_length_aa': median_len if median_len is not None else '',
        'n_distinct_cluster07': len(clusters), 'n_clusters_known_of': n_clusters_known,
        'pct_distinct_clusters': f'{pct_distinct:.3f}' if pct_distinct is not None else '',
        'n_phaC_contigs': n_phac_contigs, 'max_phaC_on_single_contig': max_on_one_contig,
        'n_phaC_contigs_known_of': n_scaffold_known,
        'checkm_completeness': comp, 'checkm_contamination': contam,
        'legitimacy_call': call,
    })

rows_out.sort(key=lambda r: -r['n_targets_corrected'])
out_path = OUT / 'phac_multicopy_legitimacy_audit.tsv'
with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)
print(f'\nwrote {out_path} ({len(rows_out):,} genomes)')

# ---------------------------------------------------------------------
# 6. summary
# ---------------------------------------------------------------------
call_counts = Counter(r['legitimacy_call'] for r in rows_out)
print(f'\nLegitimacy call breakdown (n={len(rows_out):,} genomes with n_phaC>={MIN_COPIES}):')
for call, n in call_counts.most_common():
    print(f'  {call:20s} {n:5,d}  ({100*n/len(rows_out):.1f}%)')

print(f'\nBy assembly type:')
for is_mag, label in [(True, 'MAG'), (False, 'isolate')]:
    subset = [r for r in rows_out if r['is_mag'] == is_mag]
    if not subset:
        continue
    cc = Counter(r['legitimacy_call'] for r in subset)
    print(f'  {label} (n={len(subset):,}): ' + ', '.join(f'{k}={v} ({100*v/len(subset):.0f}%)' for k, v in cc.most_common()))

print(f'\nBy copy-count bucket:')
buckets = [(3, 4, '3-4'), (5, 9, '5-9'), (10, 999, '10+')]
for lo, hi, label in buckets:
    subset = [r for r in rows_out if lo <= r['n_targets_corrected'] <= hi]
    if not subset:
        continue
    cc = Counter(r['legitimacy_call'] for r in subset)
    print(f'  {label:6s} copies (n={len(subset):,}): ' + ', '.join(f'{k}={v} ({100*v/len(subset):.0f}%)' for k, v in cc.most_common()))

print(f'\nTop 10 studies by n likely_legit genomes:')
legit_studies = Counter(r['study_id'] for r in rows_out if r['legitimacy_call'] == 'likely_legit')
for st, n in legit_studies.most_common(10):
    print(f'  {st:15s} {n}')
