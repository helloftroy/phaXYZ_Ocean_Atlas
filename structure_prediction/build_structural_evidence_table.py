"""Collapses the raw, multi-hit-per-query output of run_foldseek_search.sh
into one best-hit row per candidate, for the candidate -> bona fide PhaC
structural-evidence question in PHA_CLEAN_RESULTS.md section 6.

Two things the raw foldseek TSV does NOT give you directly, both handled here:
  1. Multiple hits per candidate (one per reference structure that cleared
     -e 10) -- this picks the single best one per candidate_id.
  2. Candidates with NO hit at all, even at that loose -e 10 -- foldseek's
     TSV just omits these rows entirely, which would silently drop them
     from any downstream distribution/count. Cross-checked against
     fold_manifest.tsv (the full candidate list from build_fold_inputs.py)
     so every candidate gets a row, status='no_hit' included.

Best-hit ranking: highest alntmscore (TM-score normalized by the
structural alignment length itself, not by query or target length
separately) wins, ties broken by bits. alntmscore was picked as the single
ranking criterion because it does not favor either long queries (like
qtmscore can, by rewarding coverage of a long query) or long references
(like ttmscore can) -- it is the most "even-handed" single number for
picking the best-matching reference structure. qtmscore/ttmscore are still
kept as separate reported columns since they answer a different, also
useful question (how much of THIS candidate's own structure looks like a
PhaC, vs. how much of the reference it matched looks recovered) -- do not
assume alntmscore-best is also qtmscore-best or ttmscore-best for a given
candidate; it usually is not exactly the same hit.

Usage (after run_foldseek_search.sh has produced its output):
    python structure_prediction/build_structural_evidence_table.py

QUERY_SETS is discovered from fold_manifest.tsv itself (every distinct
'set' value except 'reference', which is the Foldseek DB side, not a
query set) rather than hardcoded -- confirmed live (2026-09-26) that a
hardcoded list is a real, recurring failure mode: this script silently
produced nothing for all_phac_dedup even after the corresponding
ESMFold+Foldseek run had completed, because that literal name was never
added to the list here (fold_manifest.tsv itself had the same gap --
see add_all_phac_dedup_to_manifest.py). Auto-discovery means adding a
new fold set in the future only requires registering it in
fold_manifest.tsv, not also remembering to edit this file.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FOLDSEEK_OUT = ROOT / 'foldseek_out'
MANIFEST = ROOT / 'fold_manifest.tsv'
OUT_PATH = FOLDSEEK_OUT / 'structural_evidence_best_hit.tsv'


def discover_query_sets():
    sets = []
    seen = set()
    with open(MANIFEST) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            s = row['set']
            if s != 'reference' and s not in seen:
                seen.add(s)
                sets.append(s)
    return sets


QUERY_SETS = discover_query_sets()

OUT_COLUMNS = [
    'candidate_id', 'group', 'status', 'best_reference',
    'qtmscore', 'ttmscore', 'alntmscore', 'aligned_length',
    'qcov_struct', 'tcov_struct', 'structural_pident',
    'evalue', 'bits', 'prob',
]


def load_manifest_candidates():
    """All candidate IDs per group, from fold_manifest.tsv -- the source of
    truth for who SHOULD have a row, independent of whether foldseek found
    a hit for them."""
    by_group = {s: [] for s in QUERY_SETS}
    with open(MANIFEST) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            if row['set'] in by_group:
                by_group[row['set']].append(row['target_id'])
    return by_group


def load_best_hits(query_set):
    """query_set's raw foldseek TSV -> {candidate_id: best hit dict}."""
    path = FOLDSEEK_OUT / f'{query_set}_vs_reference.tsv'
    if not path.exists():
        return None  # caller distinguishes "not run yet" from "ran, zero hits"
    fieldnames = ['query', 'target', 'evalue', 'bits', 'prob', 'alntmscore',
                  'qtmscore', 'ttmscore', 'lddt', 'pident', 'qcov', 'tcov',
                  'qlen', 'tlen', 'alnlen']
    best = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t', fieldnames=fieldnames)
        for row in r:
            cid = row['query']
            score = float(row['alntmscore'])
            if cid not in best or score > float(best[cid]['alntmscore']):
                best[cid] = row
    return best


def main():
    by_group = load_manifest_candidates()
    rows = []
    for group in QUERY_SETS:
        candidates = by_group.get(group, [])
        best_hits = load_best_hits(group)
        if best_hits is None:
            print(f'WARNING: {FOLDSEEK_OUT / (group + "_vs_reference.tsv")} not found -- '
                  f'run structure_prediction/run_foldseek_search.sh first (for {group}). Skipping this group.')
            continue
        n_hit = n_no_hit = 0
        for cid in candidates:
            hit = best_hits.get(cid)
            if hit is None:
                n_no_hit += 1
                rows.append({
                    'candidate_id': cid, 'group': group, 'status': 'no_hit',
                    'best_reference': '', 'qtmscore': '', 'ttmscore': '',
                    'alntmscore': '', 'aligned_length': '', 'qcov_struct': '',
                    'tcov_struct': '', 'structural_pident': '',
                    'evalue': '', 'bits': '', 'prob': '',
                })
                continue
            n_hit += 1
            rows.append({
                'candidate_id': cid, 'group': group, 'status': 'hit',
                'best_reference': hit['target'],
                'qtmscore': hit['qtmscore'], 'ttmscore': hit['ttmscore'],
                'alntmscore': hit['alntmscore'], 'aligned_length': hit['alnlen'],
                'qcov_struct': hit['qcov'], 'tcov_struct': hit['tcov'],
                'structural_pident': hit['pident'],
                'evalue': hit['evalue'], 'bits': hit['bits'], 'prob': hit['prob'],
            })
        print(f'{group}: {n_hit} with a best hit, {n_no_hit} with no hit even at -e 10 '
              f'(of {len(candidates)} total per fold_manifest.tsv)')

    if not rows:
        print('Nothing written -- no foldseek output found yet for any query set.')
        return

    FOLDSEEK_OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=OUT_COLUMNS, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {OUT_PATH} ({len(rows)} rows)')


if __name__ == '__main__':
    main()
