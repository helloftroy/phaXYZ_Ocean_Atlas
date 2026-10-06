"""Read out the structural verdict on the 255 missed-synthase finds.

Run after both cluster jobs have finished:
    sbatch cluster/run_novel_phac_structure_prediction.sbatch
    sbatch --export=ALL,QUERY_SETS=novel_phac_candidates cluster/run_foldseek.sbatch

Sequence evidence already places these in the phaC family (section 13.3);
structure is here to answer two things it cannot.

  Is the fold real?  Foldseek alntmscore against the folded audited
    reference set, plus ESMFold's own mean pLDDT. A high TM-score to a
    validated synthase is fold-level evidence independent of every
    sequence check already run.
  Is the active site real, and is the protein whole?  The geometric triad
    search (catalytic_domain/find_structural_triad.py, same <=4.5A
    two-leg cutoff as the divergent35 audit) says whether Cys, Asp and
    His actually meet in space. tcov_struct -- how much of the matched
    REFERENCE the candidate covers -- is the truncation test: a complete
    compact synthase aligns across most of its reference, while a
    truncated gene call lights up a high TM-score over a fraction of it.

Nothing here re-decides family membership; it sorts the set into
"complete compact synthase", "truncated gene call", and "neither".

Output: structure_prediction/novel_phac_structural_verdicts.tsv

Usage:
    python structure_prediction/analyze_novel_phac_structures.py
"""
import csv
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'catalytic_domain'))
from find_structural_triad import TIGHT_CUTOFF_A, find_triad, parse_pdb_residues  # noqa: E402

HERE = ROOT / 'structure_prediction'
FASTA = HERE / 'novel_phac_candidates.faa'
PDB_DIR = HERE / 'esmfold_out/novel_phac_candidates'
FOLDSEEK = HERE / 'foldseek_out/novel_phac_candidates_vs_reference.tsv'
ANNOTATIONS = ROOT / 'phac_recovery/novel_phac_new_candidates.tsv'
OUT = HERE / 'novel_phac_structural_verdicts.tsv'

# A candidate counts as structurally confirmed on fold when its best reference
# hit clears this TM-score. 0.5 is the conventional "same fold" threshold and is
# the one every other structural call in this project uses.
TM_SAME_FOLD = 0.5
# Below this share of the matched reference covered, the candidate is aligning to
# part of a synthase rather than all of one -- the truncated-gene-call signature.
TCOV_COMPLETE = 0.75

def as_percent(plddt):
    """pLDDT reaches this project on both scales -- the PDB_files/ structures
    store it 0-1, ESMFold's own B-factor column is 0-100. Report one scale."""
    return plddt * 100 if plddt <= 1.0 else plddt


if not PDB_DIR.is_dir() or not any(PDB_DIR.glob('*.pdb')):
    raise SystemExit(f'No PDB files in {PDB_DIR} -- run cluster/run_novel_phac_structure_prediction.sbatch first.')

sequences, name, buf = {}, None, []
with open(FASTA) as f:
    for line in f:
        if line.startswith('>'):
            if name:
                sequences[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        sequences[name] = ''.join(buf)

annotations = {}
with open(ANNOTATIONS, newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        annotations[row['gene_id']] = row

# Best Foldseek hit per candidate, ranked by alntmscore then bits -- same rule as
# build_structural_evidence_table.py, so the two tables stay comparable.
best_hit = {}
if FOLDSEEK.exists():
    fields = ['query', 'target', 'evalue', 'bits', 'prob', 'alntmscore', 'qtmscore', 'ttmscore',
              'lddt', 'pident', 'qcov', 'tcov', 'qlen', 'tlen', 'alnlen']
    with open(FOLDSEEK, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t', fieldnames=fields):
            gene = row['query'].removesuffix('.pdb')
            prev = best_hit.get(gene)
            if prev is None or (float(row['alntmscore']), float(row['bits'])) > \
                    (float(prev['alntmscore']), float(prev['bits'])):
                best_hit[gene] = row
    print(f'{len(best_hit):,} candidates have a Foldseek hit')
else:
    print(f'{FOLDSEEK} not found -- fold evidence will be blank. Run cluster/run_foldseek.sbatch '
          f'with QUERY_SETS=novel_phac_candidates to fill it in.')

FIELDS = ['gene_id', 'genome', 'tier', 'pident', 'gene_len', 'gtdb_phylum', 'gtdb_genus', 'ecosystem_type',
          'mean_plddt', 'nucleophile', 'triad_tight', 'nuc_his_dist_A', 'his_asp_dist_A',
          'triad_min_plddt', 'best_reference', 'alntmscore', 'qtmscore', 'lddt', 'qcov_struct', 'tcov_struct',
          'verdict']
rows = []
for gene_id, seq in sorted(sequences.items()):
    pdb = PDB_DIR / f'{gene_id}.pdb'
    if not pdb.exists():
        continue
    residues = parse_pdb_residues(pdb)
    triad = find_triad(seq, residues, nucleophile='C')
    nucleophile = 'C'
    tight = triad is not None and triad['nuc_his_dist_A'] <= TIGHT_CUTOFF_A and triad['his_asp_dist_A'] <= TIGHT_CUTOFF_A
    if not tight:
        for alt_nuc in ('S', 'T'):
            alt = find_triad(seq, residues, nucleophile=alt_nuc)
            if alt and alt['nuc_his_dist_A'] <= TIGHT_CUTOFF_A and alt['his_asp_dist_A'] <= TIGHT_CUTOFF_A:
                triad, nucleophile, tight = alt, alt_nuc, True
                break

    hit = best_hit.get(gene_id)
    tm = float(hit['alntmscore']) if hit else None
    tcov = float(hit['tcov']) if hit else None
    if tm is None:
        verdict = 'no fold evidence'
    elif tm < TM_SAME_FOLD:
        verdict = 'fold does not match a synthase'
    elif not tight:
        verdict = 'synthase fold, no closed active site'
    elif tcov is not None and tcov < TCOV_COMPLETE:
        verdict = 'synthase fold, partial -- likely truncated gene call'
    else:
        verdict = 'complete synthase, closed active site'

    annotation = annotations.get(gene_id, {})
    rows.append({
        'gene_id': gene_id, 'genome': annotation.get('genome', ''), 'tier': annotation.get('tier', ''),
        'pident': annotation.get('pident', ''), 'gene_len': len(seq),
        'gtdb_phylum': annotation.get('gtdb_phylum', ''), 'gtdb_genus': annotation.get('gtdb_genus', ''),
        'ecosystem_type': annotation.get('ecosystem_type', ''),
        'mean_plddt': f"{as_percent(statistics.mean(r['plddt'] for r in residues.values())):.1f}" if residues else '',
        'nucleophile': nucleophile if tight else '',
        'triad_tight': str(tight),
        'nuc_his_dist_A': f"{triad['nuc_his_dist_A']:.2f}" if triad else '',
        'his_asp_dist_A': f"{triad['his_asp_dist_A']:.2f}" if triad else '',
        'triad_min_plddt': f"{as_percent(min(triad['nuc_plddt'], triad['asp_plddt'], triad['his_plddt'])):.1f}" if triad else '',
        'best_reference': hit['target'].removesuffix('.pdb') if hit else '',
        'alntmscore': hit['alntmscore'] if hit else '', 'qtmscore': hit['qtmscore'] if hit else '',
        'lddt': hit['lddt'] if hit else '', 'qcov_struct': hit['qcov'] if hit else '',
        'tcov_struct': hit['tcov'] if hit else '', 'verdict': verdict,
    })

with open(OUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=FIELDS)
    writer.writeheader()
    writer.writerows(rows)
print(f'wrote {OUT} -- {len(rows)} of {len(sequences)} candidates folded')

print('\nverdicts:')
for verdict, n in Counter(r['verdict'] for r in rows).most_common():
    print(f'  {n:4d}  {verdict}')
confirmed = [r for r in rows if r['verdict'] == 'complete synthase, closed active site']
if confirmed:
    print(f'\nconfirmed set: {len(confirmed)} genes in {len({r["genome"] for r in confirmed})} genomes')
    print('  phylum:', Counter(r['gtdb_phylum'] or 'unassigned' for r in confirmed).most_common(5))
    print('  habitat:', Counter(r['ecosystem_type'] or 'unknown' for r in confirmed).most_common(5))
    print('  tier:', dict(Counter(r['tier'] for r in confirmed)))
    print(f'  mean pLDDT: median {statistics.median(float(r["mean_plddt"]) for r in confirmed):.1f}')
    print(f'  TM-score to nearest reference: median '
          f'{statistics.median(float(r["alntmscore"]) for r in confirmed):.2f}')
