"""Builds the mmseqs query set for the novel-phaC recovery search
(searching the 7,826 qualifying_genomes.txt genomes -- no phaC hit, but
>5 other PHA-pathway families present -- for a divergent synthase our
original sequence search missed): triad-complete, QC-passing phaC
targets, one representative sequence per phaC_cluster0.7 (70%-identity)
cluster -- avoids querying with many near-redundant copies of the same
well-established paralog while still covering every distinct validated
phaC lineage this project has confirmed.

Query universe is every member of phaC_cluster0.7_cluster.tsv (the same
~128k-target population every other cluster0.7-based figure in this
project draws from -- see figures/scripts/_triad_filter.py and its
callers), not phaC_unique_targets_with_metadata.tsv's larger 191,793
raw-recruited-target count -- the gap is this project's established
downstream QC/length-floor filtering, already baked into which targets
made it into cluster0.7 in the first place.

Triad-completeness: the same hybrid check figures/scripts/_triad_filter.py
uses for the circos diagrams (structural geometry via
catalytic_domain/find_structural_triad.py where a folded PDB exists --
2,108 available locally -- alignment-column catalytic_domain/
phac_catalytic_triad.tsv fallback otherwise) -- reused directly, not
reimplemented, so this query set is defined the same way "triad-complete"
means everywhere else in this project.

Cluster representative choice: for each cluster0.7 group, prefer the
mmseqs-assigned representative id itself if it survives the triad+QC
filter (keeps the query set aligned with the id every other cluster0.7
figure already treats as "the" representative); otherwise fall back to
the lexicographically smallest surviving member of that cluster
(deterministic, not biologically meaningful -- any surviving member
represents the cluster equally well for a presence/absence recovery
search).

Usage:
  .venv/bin/python phac_recovery/build_novel_phac_query_set.py
"""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'figures' / 'scripts'))
import _phac_qc
from _triad_filter import load_triad_complete_set

FA = ROOT / 'PHA_bioprospecting' / 'omdb_search' / 'results'
CLUSTER_TSV = FA / 'phaC_cluster0.7_cluster.tsv'
SEQ_FAA = FA / 'phaC_cluster_sequences.faa'
OUT_FAA = Path(__file__).resolve().parent / 'novel_phac_query_representatives.faa'

# 1. full cluster0.7 population + membership (member -> its cluster's representative id)
member_to_rep = {}
with open(CLUSTER_TSV, newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) < 2:
            continue
        rep, member = row[0], row[1]
        member_to_rep[member] = rep

all_target_ids = set(member_to_rep)
print(f'{len(all_target_ids)} total phaC targets in phaC_cluster0.7_cluster.tsv')

# 2. QC filter (same exclusion list every figure-generating script uses)
bad = _phac_qc.load_bad_targets()
qc_pass = {t for t in all_target_ids if t not in bad}
print(f'{len(qc_pass)} pass QC ({len(all_target_ids) - len(qc_pass)} excluded by load_bad_targets)')

# 3. triad-complete filter (hybrid structural + alignment-column, same as the circos figures)
triad_complete = load_triad_complete_set(qc_pass)
print(f'{len(triad_complete)} are triad-complete (QC-passing)')

# 4. group by cluster0.7, pick one representative per cluster
clusters: dict[str, list[str]] = {}
for t in triad_complete:
    clusters.setdefault(member_to_rep[t], []).append(t)
print(f'{len(clusters)} distinct phaC_cluster0.7 clusters represented among triad-complete targets')

chosen: dict[str, str] = {}
n_fallback = 0
for rep, members in clusters.items():
    if rep in members:
        chosen[rep] = rep
    else:
        chosen[rep] = sorted(members)[0]
        n_fallback += 1
print(f'{n_fallback}/{len(clusters)} clusters used a fallback representative '
      f'(the mmseqs-assigned rep itself was not triad-complete/QC-passing)')

query_ids = set(chosen.values())
print(f'{len(query_ids)} final query sequences')

# 5. pull sequences for the chosen representatives
seqs: dict[str, str] = {}
cur_id, cur_seq = None, []
with open(SEQ_FAA) as f:
    for line in f:
        if line.startswith('>'):
            if cur_id in query_ids:
                seqs[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id in query_ids:
        seqs[cur_id] = ''.join(cur_seq)

missing = query_ids - set(seqs)
if missing:
    sys.exit(f'{len(missing)} query target_ids missing from {SEQ_FAA} -- unexpected, aborting.')

OUT_FAA.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_FAA, 'w') as f:
    for t in sorted(seqs):
        f.write(f'>{t}\n')
        seq = seqs[t]
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + '\n')

print(f'Wrote {len(seqs)} query sequences -> {OUT_FAA}')
