"""Registers all_phac_dedup_representatives.faa's 15,880 sequences into
fold_manifest.tsv, set='all_phac_dedup'.

build_all_phac_dedup_set.py (section 9.8) built that FASTA as a fourth,
independent, comprehensive structural-atlas input set -- but never wrote
its own IDs into fold_manifest.tsv the way build_fold_inputs.py did for
the original three (uncertain/positive_control/reference). Confirmed
live (2026-09-26): fold_manifest.tsv had zero all_phac_dedup rows even
after the ESMFold array job for that set had been run on the cluster --
build_structural_evidence_table.py's load_manifest_candidates() reads
fold_manifest.tsv as its source of truth for "who SHOULD have a row", so
with no rows registered it silently produced nothing for this set
(build_structural_evidence_table.py's own QUERY_SETS list also needed a
matching fix -- see that script).

Idempotent: only appends target_ids not already present under ANY set
(668 of the 15,880 already appear as uncertain/positive_control/
reference representatives -- kept as their own separate all_phac_dedup
row too, consistent with each set being folded/searched independently
in its own esmfold_out/<set>/ directory; not deduped against the other
three).

Usage:
    python structure_prediction/add_all_phac_dedup_to_manifest.py
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / 'fold_manifest.tsv'
FASTA = ROOT / 'all_phac_dedup_representatives.faa'
SET_NAME = 'all_phac_dedup'

existing_for_set = set()
all_lines = []
if MANIFEST.exists():
    with open(MANIFEST) as f:
        header = f.readline()
        all_lines.append(header)
        for line in f:
            all_lines.append(line)
            tid, s = line.rstrip('\n').split('\t')
            if s == SET_NAME:
                existing_for_set.add(tid)
else:
    all_lines.append('target_id\tset\n')

new_ids = []
with open(FASTA) as f:
    for line in f:
        if line.startswith('>'):
            tid = line[1:].split()[0]
            if tid not in existing_for_set:
                new_ids.append(tid)

print(f'{len(existing_for_set)} {SET_NAME} rows already in {MANIFEST}, {len(new_ids)} new to add')

if new_ids:
    with open(MANIFEST, 'w') as out:
        out.writelines(all_lines)
        for tid in new_ids:
            out.write(f'{tid}\t{SET_NAME}\n')
    print(f'wrote {MANIFEST}')
else:
    print('nothing to do')
