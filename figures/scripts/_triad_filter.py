"""Hybrid triad-completeness check, used to filter the circos diagrams
down to only triad-complete phaC copies (per direct request, on top of
the existing qtmscore>=0.5 structural-confidence filter): the structural
geometric check (catalytic_domain/find_structural_triad.py, tight <=4.5A
cutoff, canonical Cys or alternative Ser/Thr nucleophile) where a folded
PDB exists, falling back to the older alignment-column triad_complete
flag (catalytic_domain/phac_catalytic_triad.tsv) where it doesn't.

The structural check takes precedence, rather than using the alignment-
column flag alone for consistency, because sections 9.11/9.12 already
found the alignment-column method wrong in real cases -- e.g. CARD22-1's
708aa phaC/phaZ-ambiguous copy: the alignment-column method marks it
triad_complete=False (hmmalign's glocal Viterbi algorithm tends to
collapse divergent sequences into insert states, missing real match-
column residues), but the structure itself shows a textbook-tight
Cys-Asp-His triad (3.19A/2.73A). The fallback exists because most
circos-plotted targets do NOT have a folded PDB (structural_evidence_
best_hit.tsv only covers the uncertain/positive_control/no_hmm_triad_
support sets, not every genome's own paralogs) -- the alignment-column
flag has full dataset coverage, unlike structure.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / 'catalytic_domain'))
from find_structural_triad import structural_triad_complete_batch

ROOT = Path(__file__).resolve().parent.parent.parent


def load_triad_complete_set(target_ids):
    """Returns the subset of target_ids that are triad-complete by the
    hybrid check above."""
    target_ids = set(target_ids)
    structural = structural_triad_complete_batch(target_ids)

    alignment_triad = {}
    with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['target_id'] in target_ids:
                alignment_triad[row['target_id']] = row['triad_complete'] == 'True'

    good = set()
    n_fallback = 0
    for t in target_ids:
        if t in structural:
            if structural[t]:
                good.add(t)
        else:
            n_fallback += 1
            if alignment_triad.get(t, False):
                good.add(t)
    print(f'triad-complete filter: {len(structural)} targets resolved structurally, '
          f'{n_fallback} fell back to the alignment-column flag; {len(good)}/{len(target_ids)} kept')
    return good
