"""Builds a structure-prediction input covering essentially all verified
phaC diversity (not a hand-picked subset like uncertain_cluster_
representatives.faa/build_uncertain_set_v2.py) -- the goal is a reference
structural atlas broad enough to support future figures (e.g. the
sequence-similarity-vs-structural-similarity plot planned next), not to
resolve any one specific question.

Two cuts, both requested directly rather than invented here:
  1. Drop short sequences (<MIN_LEN aa) -- not structurally interpretable
     regardless of source, same floor used throughout this project.
  2. Drop near-duplicates (>=REDUNDANCY_IDENTITY% identity to another
     candidate) -- folding two 99%-identical sequences wastes GPU time for
     no new structural information. No length ceiling is applied here
     (unlike build_uncertain_set_v2.py's two dataset-wide sources) --
     this run is meant to be comprehensive, not filtered to "typical"
     lengths.

No mmseqs2/CD-HIT binary is available in this environment, so the
redundancy cut is built directly from data this project already has:
catalytic_domain/selfsearch/hits.tsv, the same all-vs-all phaC self-search
used throughout section 9 (re-filtered against the current bad-target
list, same discipline as build_selfsearch_derived_tables.py). Two targets
are called near-duplicates if pident>=REDUNDANCY_IDENTITY AND both
qcov>=REDUNDANCY_COV and tcov>=REDUNDANCY_COV (near-full-length match on
both sides, not just a locally-identical fragment -- the same two-sided
coverage discipline mmseqs clustering itself uses, just implemented here
as a union-find over the self-search graph instead of a fresh mmseqs run).
Connected components of this graph are collapsed to one representative
each (longest sequence in the component, deterministic).

IMPORTANT CAVEAT on this redundancy cut: the self-search was run with
mmseqs2 --max-seqs 20 (see build_selfsearch_derived_tables.py's docstring)
-- each query only has its top ~20 hits recorded, not a truly exhaustive
all-vs-all. For most candidates this does not matter (a 95%+-identical
duplicate is almost always in the query's top 20 hits by construction),
but a candidate whose top 20 hits are all UNRELATED sequences and whose
one true near-duplicate happens to rank 21st would not get merged here.
This makes the redundancy cut a slight UNDER-estimate (a few extra
near-duplicate pairs may survive as separate entries) rather than an
over-estimate that could wrongly collapse genuinely distinct sequences --
the safer direction for a structural atlas meant to be comprehensive.

Usage:
    python structure_prediction/build_all_phac_dedup_set.py

Outputs:
    structure_prediction/all_phac_dedup_representatives.faa
    structure_prediction/all_phac_dedup_clusters.tsv (representative_id -> all collapsed member target_ids)
"""
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = Path(__file__).resolve().parent

sys.path.insert(0, str(ROOT / 'figures/scripts'))
import _phac_qc

MIN_LEN = 150
REDUNDANCY_IDENTITY = 95.0
REDUNDANCY_COV = 0.90

bad_targets = _phac_qc.load_bad_targets()

# ---------------------------------------------------------------------
# 1. full verified candidate universe (corrected join, current exclusion list)
# ---------------------------------------------------------------------
genome_targets_raw = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
all_candidates = set()
for g, ts in genome_targets_raw.items():
    all_candidates |= (ts - bad_targets)
print(f'{len(all_candidates):,} verified candidate target_ids (corrected)')

target_lengths = pickle.load(open('/tmp/all_phac_target_lengths.pkl', 'rb'))
length_filtered = {t for t in all_candidates if target_lengths.get(t, 0) >= MIN_LEN}
print(f'{len(length_filtered):,} after dropping <{MIN_LEN}aa '
      f'({len(all_candidates) - len(length_filtered):,} dropped)')

# ---------------------------------------------------------------------
# 2. union-find over the self-search graph, edges = near-duplicate pairs
# ---------------------------------------------------------------------
parent = {t: t for t in length_filtered}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb


n_edges = 0
with open(ROOT / 'catalytic_domain/selfsearch/hits.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        q, t = row['query'], row['target']
        if q == t or q not in length_filtered or t not in length_filtered:
            continue
        if float(row['pident']) >= REDUNDANCY_IDENTITY and float(row['qcov']) >= REDUNDANCY_COV and float(row['tcov']) >= REDUNDANCY_COV:
            union(q, t)
            n_edges += 1
print(f'{n_edges:,} near-duplicate edges (>={REDUNDANCY_IDENTITY:.0f}% identity, >={REDUNDANCY_COV:.0%} coverage both sides)')

components = defaultdict(list)
for t in length_filtered:
    components[find(t)].append(t)

n_components = len(components)
n_multi = sum(1 for members in components.values() if len(members) > 1)
print(f'{n_components:,} redundancy clusters ({n_multi:,} with >1 member, '
      f'{len(length_filtered) - n_components:,} sequences collapsed away)')

representatives = {}
for root, members in components.items():
    rep = max(members, key=lambda t: target_lengths.get(t, 0))
    representatives[rep] = sorted(members)

print(f'\n{len(representatives):,} final representative sequences to fold')

# ---------------------------------------------------------------------
# 3. write cluster mapping + FASTA
# ---------------------------------------------------------------------
cluster_out = OUT / 'all_phac_dedup_clusters.tsv'
with open(cluster_out, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['representative_id', 'n_members', 'member_target_ids'])
    for rep in sorted(representatives):
        members = representatives[rep]
        w.writerow([rep, len(members), ','.join(members)])
print(f'wrote {cluster_out}')

sequences = {}
cur_id, cur_seq = None, []
with open(FA / 'phaC_cluster_sequences.faa') as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in representatives:
                sequences[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in representatives:
        sequences[cur_id] = ''.join(cur_seq)

missing = set(representatives) - set(sequences)
if missing:
    print(f'WARNING: {len(missing)} representative target_ids not found in phaC_cluster_sequences.faa')

fasta_out = OUT / 'all_phac_dedup_representatives.faa'
n_written = 0
with open(fasta_out, 'w') as out:
    for rep in sorted(representatives):
        seq = sequences.get(rep)
        if not seq:
            continue
        out.write(f'>{rep}\n')
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + '\n')
        n_written += 1
print(f'wrote {fasta_out} ({n_written} sequences)')

lens = sorted(target_lengths.get(t, 0) for t in representatives if t in sequences)
if lens:
    print(f'\nlength range: {lens[0]}-{lens[-1]}aa, median {lens[len(lens)//2]}aa')
