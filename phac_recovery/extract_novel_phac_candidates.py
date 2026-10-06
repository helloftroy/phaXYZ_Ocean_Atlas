"""Tier the mmseqs hits from the missed-synthase search and extract the
candidate protein sequences.

Input:  phac_recovery/novel_phac_mmseqs_hits.tsv  (cluster/run_novel_phac_mmseqs_search.sbatch)
        phac_recovery/novel_phac_target_proteomes.faa  (only needed for the FASTA output)
Output: phac_recovery/novel_phac_candidates.tsv   one row per candidate gene, best hit + tier
        phac_recovery/novel_phac_candidates.faa   sequences, header = gene id + tier

Tier is set by each gene's best hit (highest bit score) to a validated phaC:
  A  >=50% identity, >=80% of both query and gene aligned
  B  30-50% identity, >=80% of both
  C  <30% identity, >=80% of both
  D  >=80% of the gene aligned but not of the query (short gene)
  E  partial / domain-only

Standard library only. Without the proteome file it writes the table and skips the FASTA.

Usage:
    python phac_recovery/extract_novel_phac_candidates.py
"""
import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
HITS = HERE / 'novel_phac_mmseqs_hits.tsv'
PROTEOMES = HERE / 'novel_phac_target_proteomes.faa'
OUT_TSV = HERE / 'novel_phac_candidates.tsv'
OUT_FAA = HERE / 'novel_phac_candidates.faa'


def tier_of(pid, qcov, tcov):
    if qcov >= 0.8 and tcov >= 0.8:
        return 'A' if pid >= 50 else 'B' if pid >= 30 else 'C'
    return 'D' if tcov >= 0.8 else 'E'


best = {}
n_queries = Counter()
with open(HITS) as f:
    # columns: query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen
    for r in csv.reader(f, delimiter='\t'):
        gene, bits = r[1], float(r[11])
        n_queries[gene] += 1
        if gene not in best or bits > best[gene]['bits']:
            qlen, tlen = int(r[12]), int(r[13])
            best[gene] = dict(
                best_query=r[0], pident=float(r[2]), evalue=r[10], bits=bits, qlen=qlen, tlen=tlen,
                qcov=(int(r[7]) - int(r[6]) + 1) / qlen, tcov=(int(r[9]) - int(r[8]) + 1) / tlen)

for gene, b in best.items():
    b['tier'] = tier_of(b['pident'], b['qcov'], b['tcov'])

with open(OUT_TSV, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene_id', 'genome', 'tier', 'best_query', 'pident', 'qcov', 'tcov', 'query_len', 'gene_len',
                'evalue', 'bits', 'n_queries_hit'])
    for gene in sorted(best, key=lambda g: (best[g]['tier'], -best[g]['bits'])):
        b = best[gene]
        w.writerow([gene, gene.split('-scaffold_')[0], b['tier'], b['best_query'], f"{b['pident']:.1f}",
                    f"{b['qcov']:.3f}", f"{b['tcov']:.3f}", b['qlen'], b['tlen'], b['evalue'], f"{b['bits']:.0f}",
                    n_queries[gene]])
counts = Counter(b['tier'] for b in best.values())
print(f'{len(best):,} candidate genes -> {OUT_TSV.name}: ' + ', '.join(f'{t}={counts[t]:,}' for t in 'ABCDE'))

if not PROTEOMES.exists():
    print(f'{PROTEOMES.name} not found here, so no FASTA written (run this on the cluster for the sequences).')
else:
    n_written = 0
    keep = False
    with open(PROTEOMES) as fin, open(OUT_FAA, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                gene = line[1:].split()[0]
                keep = gene in best
                if keep:
                    fout.write(f">{gene} tier={best[gene]['tier']}\n")
                    n_written += 1
            elif keep:
                fout.write(line)
    print(f'{n_written:,} of {len(best):,} sequences -> {OUT_FAA.name}')
