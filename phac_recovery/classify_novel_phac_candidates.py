"""Decide which of the triad-complete missed-synthase candidates are actually
new, and describe the ones that are.

The mmseqs search (cluster/run_novel_phac_mmseqs_search.sbatch) searched the
whole proteomes of PHA-rich, apparently phaC-negative genomes with validated
triad-complete phaC queries, and check_novel_phac_candidate_triads.py scored
every hit for the catalytic triad. Most triad-complete hits are NOT new: they
are proteins the atlas already holds, surfaced again because the genome was
only ever recorded phaC-negative for a bookkeeping reason. This script
separates the three cases by matching each candidate's sequence against the
phaC target sequences (NR100, so amino-acid-exact; the proteome FASTA writes a
terminal '*' where the target FASTA writes 'X', hence the normalisation):

  known_qc_passing   exactly a phaC target that PASSED QC. The genome really is
                     phaC-positive and was mis-recorded by the 5-genome
                     metadata cap (see PHA_CLEAN_RESULTS.md section 2).
  known_qc_removed   exactly a phaC target that QC DROPPED, because the single
                     reference it best matched is on the bad-reference list.
                     A QC false negative, not a discovery.
  new                not in the phaC target set at all, at any identity. These
                     are the genuine finds: proteins OMDB never assigned to the
                     phaC family, so no reference-based pass could have seen
                     them.

A 'new' gene in a genome that separately holds a QC-passing phaC is a new
paralog in an already-positive genome; it is reported but kept out of the
headline set, which is new genes in genomes with no QC-passing phaC at all.

Each gene in that headline set is then scored against the three NCBIFam class
models (TIGR01838 I / TIGR01839 II / TIGR01836 III), same convention as
catalytic_domain/classify_phac_synthase_class.py, and joined to GTDB taxonomy,
study and habitat.

Outputs
  phac_recovery/novel_phac_candidate_verdicts.tsv   every triad-complete gene
  phac_recovery/novel_phac_new_candidates.faa       the headline set, for folding
  phac_recovery/novel_phac_new_candidates.tsv       the headline set, annotated

Usage:
    python phac_recovery/classify_novel_phac_candidates.py
"""
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pyhmmer

ROOT = Path(__file__).resolve().parent.parent
HERE = ROOT / 'phac_recovery'
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
HMM_DIR = HERE / 'hmm/extra'
sys.path.insert(0, str(ROOT / 'figures' / 'scripts'))
import _phac_qc

MODELS = {'TIGR01838': 'I', 'TIGR01839': 'II', 'TIGR01836': 'III'}


def norm(seq):
    """Compare sequences the two FASTAs' way: drop stop codons and a trailing X."""
    return seq.upper().replace('*', '').rstrip('X')


def read_fasta(path):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, ''.join(buf)


# ---------------------------------------------------------------- inputs
bad = _phac_qc.load_bad_targets()

seq_to_target = {}
for target_id, seq in read_fasta(FA / 'phaC_cluster_sequences.faa'):
    seq_to_target.setdefault(norm(seq), target_id)

candidate_seq = {gene: norm(seq) for gene, seq in read_fasta(HERE / 'novel_phac_candidates.faa')}

genome_has_phac = set()
genome_in_preqc = set()
with open(ROOT / 'phaC_all_genomes_from_nr100_clusters.tsv', newline='') as f:
    r = csv.reader(f, delimiter='\t')
    next(r)
    for target_id, genome in r:
        genome_in_preqc.add(genome)
        if target_id not in bad:
            genome_has_phac.add(genome)

triads = [r for r in csv.DictReader(open(HERE / 'novel_phac_candidate_triads.tsv'), delimiter='\t')
          if r['triad_complete'] == 'True']
print(f'{len(seq_to_target):,} phaC target sequences, {len(candidate_seq):,} candidates, '
      f'{len(triads):,} of them triad complete')

# ---------------------------------------------------------------- verdicts
for row in triads:
    target = seq_to_target.get(candidate_seq[row['gene_id']])
    if target is None:
        row['verdict'] = 'new'
    elif target in bad:
        row['verdict'] = 'known_qc_removed'
    else:
        row['verdict'] = 'known_qc_passing'
    row['matched_target'] = target or ''
    row['genome_has_qc_passing_phac'] = str(row['genome'] in genome_has_phac)

counts = Counter(r['verdict'] for r in triads)
by_tier = defaultdict(Counter)
for r in triads:
    by_tier[r['verdict']][r['tier']] += 1
print('\nverdicts over all triad-complete candidates:')
for v in ('known_qc_passing', 'known_qc_removed', 'new'):
    tiers = ' '.join(f'{t}={by_tier[v][t]}' for t in 'ABCDE' if by_tier[v][t])
    print(f'  {v:18s} {counts[v]:4d} genes in {len({r["genome"] for r in triads if r["verdict"] == v}):4d} genomes   {tiers}')

headline = [r for r in triads if r['verdict'] == 'new' and r['genome'] not in genome_has_phac]
paralogs = [r for r in triads if r['verdict'] == 'new' and r['genome'] in genome_has_phac]
print(f'\nheadline set: {len(headline):,} new genes in {len({r["genome"] for r in headline}):,} genomes with no '
      f'QC-passing phaC ({len(paralogs)} further new genes sit in already-positive genomes)')
print(f'  {len({r["genome"] for r in headline} & genome_in_preqc):,} of those genomes did have a phaC-family hit '
      f'before QC, but on a different protein')

# ---------------------------------------------------------------- synthase class
alphabet = pyhmmer.easel.Alphabet.amino()
wanted = {r['gene_id'] for r in headline}
seqs = []
with pyhmmer.easel.SequenceFile(str(HERE / 'novel_phac_candidates.faa'), digital=True, alphabet=alphabet) as sf:
    for s in sf:
        name = s.name.decode() if isinstance(s.name, bytes) else s.name
        if name in wanted:
            seqs.append(s)
block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)

scores, cutoff = defaultdict(dict), {}
for acc, cls in MODELS.items():
    with pyhmmer.plan7.HMMFile(str(HMM_DIR / f'{acc}.hmm')) as hf:
        hmm = hf.read()
    cutoff[cls] = hmm.cutoffs.trusted[0]
    for hits in pyhmmer.hmmsearch([hmm], block, E=10.0):
        for hit in hits:
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            scores[name][cls] = hit.score

# ---------------------------------------------------------------- metadata
taxonomy = {}
with open(FA / 'omdb_all_genome_taxonomy.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        taxonomy[row['genome']] = row
habitat = {}
with open(ROOT / 'omdb_all_genomes_with_locations.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        habitat[row['genome']] = row

FIELDS = ['gene_id', 'genome', 'tier', 'pident', 'gene_len', 'verdict', 'matched_target',
          'genome_has_qc_passing_phac', 'genome_had_preqc_phac_hit']
with open(HERE / 'novel_phac_candidate_verdicts.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=FIELDS, extrasaction='ignore')
    w.writeheader()
    for row in sorted(triads, key=lambda r: (r['verdict'], r['tier'], -float(r['pident']))):
        row['genome_had_preqc_phac_hit'] = str(row['genome'] in genome_in_preqc)
        w.writerow(row)
print(f'\nwrote {HERE / "novel_phac_candidate_verdicts.tsv"}')

OUT_FIELDS = ['gene_id', 'genome', 'tier', 'pident', 'gene_len', 'synthase_class', 'class_score', 'class_confident',
              'score_I', 'score_II', 'score_III', 'gtdb_phylum', 'gtdb_class', 'gtdb_order', 'gtdb_family',
              'gtdb_genus', 'study_id', 'ecosystem_type', 'ecosystem_name', 'latitude_degN', 'longitude_degE']
annotated = []
for row in sorted(headline, key=lambda r: (r['tier'], -float(r['pident']))):
    sc = scores.get(row['gene_id'], {})
    best = max(sc, key=sc.get) if sc else None
    tax = taxonomy.get(row['genome'], {})
    hab = habitat.get(row['genome'], {})
    annotated.append({
        'gene_id': row['gene_id'], 'genome': row['genome'], 'tier': row['tier'],
        'pident': row['pident'], 'gene_len': row['gene_len'],
        'synthase_class': best or 'none', 'class_score': f'{sc[best]:.1f}' if best else '',
        'class_confident': str(bool(best) and sc[best] >= cutoff[best]),
        'score_I': f'{sc.get("I", 0):.1f}', 'score_II': f'{sc.get("II", 0):.1f}', 'score_III': f'{sc.get("III", 0):.1f}',
        'gtdb_phylum': tax.get('gtdb_phylum', ''), 'gtdb_class': tax.get('gtdb_class', ''),
        'gtdb_order': tax.get('gtdb_order', ''), 'gtdb_family': tax.get('gtdb_family', ''),
        'gtdb_genus': tax.get('gtdb_genus', ''), 'study_id': hab.get('study_id', ''),
        'ecosystem_type': hab.get('ecosystem_type', ''), 'ecosystem_name': hab.get('ecosystem_name', ''),
        'latitude_degN': hab.get('latitude_degN', ''), 'longitude_degE': hab.get('longitude_degE', ''),
    })
with open(HERE / 'novel_phac_new_candidates.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=OUT_FIELDS)
    w.writeheader()
    w.writerows(annotated)
print(f'wrote {HERE / "novel_phac_new_candidates.tsv"}')

with open(HERE / 'novel_phac_new_candidates.faa', 'w') as f:
    for row in annotated:
        f.write(f">{row['gene_id']} tier={row['tier']} class={row['synthase_class']} pident={row['pident']}\n")
        seq = candidate_seq[row['gene_id']]
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + '\n')
print(f'wrote {HERE / "novel_phac_new_candidates.faa"} ({len(annotated)} sequences)')

# ---------------------------------------------------------------- summary
print('\nheadline set by tier:')
for tier in 'ABCDE':
    sub = [r for r in annotated if r['tier'] == tier]
    if not sub:
        continue
    pid = [float(r['pident']) for r in sub]
    ln = [int(r['gene_len']) for r in sub]
    print(f'  {tier}: n={len(sub):3d}  identity med={statistics.median(pid):.0f}% (min {min(pid):.0f}%)  '
          f'length med={statistics.median(ln):.0f} aa  confident class call={sum(r["class_confident"] == "True" for r in sub)}')
print('\nclass:', dict(Counter(r['synthase_class'] for r in annotated)))
print('phylum:', Counter(r['gtdb_phylum'] or 'unassigned' for r in annotated).most_common(8))
print('class (GTDB):', Counter(r['gtdb_class'] or 'unassigned' for r in annotated).most_common(8))
print('habitat:', Counter(r['ecosystem_type'] or 'unknown' for r in annotated).most_common(8))
print('studies:', len({r['study_id'] for r in annotated}), Counter(r['study_id'] or 'unknown' for r in annotated).most_common(6))
