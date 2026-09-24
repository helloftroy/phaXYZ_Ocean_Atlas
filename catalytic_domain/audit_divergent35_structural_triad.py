"""Runs find_structural_triad.py's geometric search (Cys, then Ser/Thr as
an alternative nucleophile if no Cys triad is found) across every
figures/phac_divergent35_validated.tsv candidate that has a folded PDB in
PDB_files/ -- these are the most sequence-divergent (<35% identity to
their nearest reference) validated phaC candidates in the dataset, each
independently confirmed by HMM match and/or PHA-pathway gene context,
not just by structure.

Built to answer a direct question raised while picking a few divergent35
structures to render for figures/PHA_CLEAN_RESULTS.md: several of the
ones with catalytic_cys_confirmed=False (no Cys anywhere in the raw
sequence) turned out, on a first spot-check, to have a real, tight,
motif-embedded SERINE-nucleophile triad -- geometrically indistinguishable
from the canonical Cys triads found in the three deep-dive genomes. This
is a different population and a different answer from figures/
structural_no_cys_alternative_triad_audit.tsv, which ran the identical
Ser/Thr check on the 15 no_hmm_triad_support/no-cysteine candidates and
found 0/15 plausible -- that population has independent-of-structure
support (an HMM match or pathway-context flag) already in hand, so
finding real alternative-nucleophile geometry in a meaningful fraction
of them here is a materially different, more interesting result than
finding none in a population with no non-structural evidence at all.

"Tight" here means both key hydrogen-bond distances (nucleophile-atom to
His, His to Asp) are under 4.5A, roughly the spread seen among this
project's known-real triads (canonical Cys ones from the CARD22-1/
Modicisalibacter/divergent35-with-confirmed-Cys structures all land at
2.4-3.4A on both legs) -- a generous cutoff, not a strict one, chosen to
separate "plausible real triad" from the clearly-nonsensical double-digit-
Angstrom "best available" combinations the no-Cys no_hmm_triad_support
population showed throughout.

Usage:
    python catalytic_domain/audit_divergent35_structural_triad.py

Output:
    figures/divergent35_structural_triad_audit.tsv
"""
import csv
import os
from pathlib import Path

from find_structural_triad import find_triad_for_target

ROOT = Path(__file__).resolve().parent.parent
TIGHT_CUTOFF_A = 4.5

pdbs = set(f[:-4] for f in os.listdir(ROOT / 'PDB_files') if f.endswith('.pdb'))

candidates = {}
with open(ROOT / 'figures/phac_divergent35_validated.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        t = row['target_id']
        if t in pdbs and t not in candidates:
            candidates[t] = row

print(f'{len(candidates)} divergent35 candidates with a folded PDB available')

out_rows = []
for t, row in sorted(candidates.items(), key=lambda kv: float(kv[1]['pident_to_nearest'])):
    best = None
    cys = find_triad_for_target(t, 'C')
    if cys is not None:
        best = ('Cys', cys)
    if best is None or best[1]['nuc_his_dist_A'] > TIGHT_CUTOFF_A or best[1]['his_asp_dist_A'] > TIGHT_CUTOFF_A:
        for nuc in ('S', 'T'):
            alt = find_triad_for_target(t, nuc)
            if alt is not None and (best is None or alt['score'] < best[1]['score']):
                best = (nuc, alt)

    if best is None:
        out_rows.append({'target_id': t, 'pident_to_nearest': row['pident_to_nearest'],
                          'confidence_tier': row['confidence_tier'], 'genus': row['genus'],
                          'catalytic_cys_confirmed_hmm': row['catalytic_cys_confirmed'],
                          'nucleophile': 'none', 'nuc_pos': '', 'asp_pos': '', 'his_pos': '',
                          'nuc_his_dist_A': '', 'his_asp_dist_A': '', 'nuc_in_motif': '', 'tight': False})
        continue
    label, d = best
    tight = d['nuc_his_dist_A'] <= TIGHT_CUTOFF_A and d['his_asp_dist_A'] <= TIGHT_CUTOFF_A
    out_rows.append({'target_id': t, 'pident_to_nearest': row['pident_to_nearest'],
                      'confidence_tier': row['confidence_tier'], 'genus': row['genus'],
                      'catalytic_cys_confirmed_hmm': row['catalytic_cys_confirmed'],
                      'nucleophile': label, 'nuc_pos': d['nuc_pos'], 'asp_pos': d['asp_pos'], 'his_pos': d['his_pos'],
                      'nuc_his_dist_A': round(d['nuc_his_dist_A'], 2), 'his_asp_dist_A': round(d['his_asp_dist_A'], 2),
                      'nuc_in_motif': d['nuc_in_motif'], 'tight': tight})

out_path = ROOT / 'figures/divergent35_structural_triad_audit.tsv'
with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter='\t')
    w.writeheader()
    w.writerows(out_rows)
print('wrote', out_path)

n_tight = sum(1 for r in out_rows if r['tight'])
print(f'\n{n_tight}/{len(out_rows)} ({100*n_tight/len(out_rows):.0f}%) have a "tight" (<= {TIGHT_CUTOFF_A}A both legs) structural triad, any nucleophile')

by_cys = {'True': [0, 0], 'False': [0, 0]}
for r in out_rows:
    k = r['catalytic_cys_confirmed_hmm']
    by_cys[k][0] += 1
    by_cys[k][1] += r['tight']
print('by HMM catalytic_cys_confirmed flag:')
for k, (n, nt) in by_cys.items():
    print(f'  confirmed={k}: {nt}/{n} tight ({100*nt/n:.0f}%)')

by_tier = {}
for r in out_rows:
    tier = r['confidence_tier'].split(':')[0]
    by_tier.setdefault(tier, [0, 0])
    by_tier[tier][0] += 1
    by_tier[tier][1] += r['tight']
print('by confidence tier:')
for tier, (n, nt) in sorted(by_tier.items(), key=lambda kv: -kv[1][1] / kv[1][0]):
    print(f'  {tier}: {nt}/{n} tight ({100*nt/n:.0f}%)')

nuc_counts = {}
for r in out_rows:
    if r['tight']:
        nuc_counts[r['nucleophile']] = nuc_counts.get(r['nucleophile'], 0) + 1
print('tight-triad nucleophile breakdown:', nuc_counts)
