"""Catalytic-triad check for the missed-synthase candidates
(phac_recovery/novel_phac_candidates.faa), using the same definition as
catalytic_domain/phac_catalytic_triad.tsv: align each sequence to
phaC_custom.hmm and read model columns 307 (Cys), 460 (Asp), 489 (His),
plus the G-x-C-x-G box at 305-309. Uses pyhmmer's hmmalign (no hmmalign
binary on this machine). Run with --validate to first re-score a sample of
already-scored validated phaC and report agreement with the existing table.

Output: phac_recovery/novel_phac_candidate_triads.tsv
"""
import csv
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyhmmer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WANTED = [305, 306, 307, 308, 309, 460, 489]
GAP = {'-', '.'}
alphabet = pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(str(HERE / 'hmm/phaC_custom.hmm')) as hf:
    HMM = hf.read()


def triad_columns(fasta, keep=None):
    seqs = []
    with pyhmmer.easel.SequenceFile(str(fasta), digital=True, alphabet=alphabet) as sf:
        for s in sf:
            name = s.name.decode() if isinstance(s.name, bytes) else s.name
            if keep is None or name in keep:
                seqs.append(s)
    out = {}
    for i in range(0, len(seqs), 2000):   # in chunks: insert columns make one giant alignment very wide
        block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs[i:i + 2000])
        msa = pyhmmer.hmmer.hmmalign(HMM, block, trim=False, all_consensus_cols=True)
        rf = msa.reference
        rf = rf.decode() if isinstance(rf, bytes) else rf
        col_of, m = {}, 0
        for j, ch in enumerate(rf):
            if ch not in '.-~':
                m += 1
                if m in WANTED:
                    col_of[m] = j
        tmsa = msa if hasattr(msa, 'alignment') and isinstance(msa.alignment[0], str) else msa.textize()
        for name, row in zip(tmsa.names, tmsa.alignment):
            name = name.decode() if isinstance(name, bytes) else name
            out[name] = {c: row[col_of[c]].upper() for c in WANTED}
    return out


def call(d):
    cys, asp, his = d[307] == 'C', d[460] == 'D', d[489] == 'H'
    return dict(cys_ok=cys, asp_ok=asp, his_ok=his, triad_complete=cys and asp and his,
                ser_thr_at_nucleophile=d[307] in 'ST',
                lipase_box_match=d[305] == 'G' and d[307] == 'C' and d[309] == 'G',
                fully_covered=all(d[c] not in GAP for c in (307, 460, 489)))


if '--validate' in sys.argv:
    known = {}
    with open(ROOT / 'catalytic_domain/phac_catalytic_triad.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            known[row['target_id']] = row['triad_complete'] == 'True'
    random.seed(0)
    sample = set(random.sample(sorted(known), 3000))
    cols = triad_columns(ROOT / 'PHA_bioprospecting/omdb_search/results/phaC_cluster_sequences.faa', sample)
    agree = Counter((known[t], call(d)['triad_complete']) for t, d in cols.items())
    n = sum(agree.values())
    print(f'validation on {n:,} already-scored proteins (existing table, this script): {dict(agree)}')
    print(f'agreement {100 * (agree[(True, True)] + agree[(False, False)]) / n:.1f}%')

cols = triad_columns(HERE / 'novel_phac_candidates.faa')
info = {r['gene_id']: r for r in csv.DictReader(open(HERE / 'novel_phac_candidates.tsv'), delimiter='\t')}
by_tier = defaultdict(Counter)
with open(HERE / 'novel_phac_candidate_triads.tsv', 'w', newline='') as f:
    fields = ['gene_id', 'genome', 'tier', 'pident', 'gene_len', 'col305', 'col306', 'col307_cys', 'col308', 'col309',
              'col460_asp', 'col489_his', 'cys_ok', 'asp_ok', 'his_ok', 'triad_complete', 'ser_thr_at_nucleophile',
              'lipase_box_match', 'fully_covered']
    w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    w.writeheader()
    for gene, d in cols.items():
        c, r = call(d), info[gene]
        w.writerow(dict(gene_id=gene, genome=r['genome'], tier=r['tier'], pident=r['pident'], gene_len=r['gene_len'],
                        col305=d[305], col306=d[306], col307_cys=d[307], col308=d[308], col309=d[309],
                        col460_asp=d[460], col489_his=d[489], **c))
        t = by_tier[r['tier']]
        t['n'] += 1
        t['triad'] += c['triad_complete']
        t['covered'] += c['fully_covered']
        t['ser_thr_full'] += c['ser_thr_at_nucleophile'] and c['asp_ok'] and c['his_ok']
        if c['triad_complete']:
            t['genomes_' + r['genome']] = 1
print(f'{len(cols):,} candidates scored -> novel_phac_candidate_triads.tsv')
for tier in sorted(by_tier):
    t = by_tier[tier]
    ng = sum(1 for k in t if k.startswith('genomes_'))
    print(f"  tier {tier}: {t['n']:6,} genes | triad complete {t['triad']:5,} ({100 * t['triad'] / t['n']:5.1f}%) in {ng:5,} genomes | "
          f"all 3 positions aligned {t['covered']:5,} | Ser/Thr + Asp + His {t['ser_thr_full']:4,}")
