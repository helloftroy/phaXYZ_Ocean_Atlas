"""Builds the three FASTA input sets for the ESMFold structure-prediction
run, per the plan in figures/PHA_CLEAN_RESULTS.md section 6.

Rather than folding all 68,424 verified phaC genomes individually, this
dedupes to one representative sequence per phaC_cluster0.7 cluster
(70%-identity groups) -- the mmseqs cluster_id IS its own representative
member's sequence ID (see phaatlas/pipeline/phylogenetics.py's docstring),
so no extra "which sequence represents this cluster" logic is needed,
just a length filter to exclude fragments (150-700aa, the same bound used
throughout this project, e.g. figures/scripts/plot_phac_divergent35_vs_depth.py).

Three sets:
  1. uncertain_cluster_representatives.faa -- one representative per
     cluster touched by a genome in any of the three "no direct evidence"
     tiers from PHA_CLEAN_RESULTS.md section 3 (>=5 other genes, 1-4
     other genes, phaC only). This is the actual target of the analysis:
     does the FOLDED STRUCTURE resemble a real PhaC catalytic domain for
     these, independent of the sequence-level evidence already tried.
  2. positive_control_representatives.faa -- a stratified sample (by
     phylum) of clusters touched ONLY by triad-complete/HMM-supported
     genomes, to confirm ESMFold + Foldseek-vs-reference actually recovers
     the right answer for genomes we're already confident about, before
     trusting it on the uncertain ones.
  3. audited_reference_set.faa -- copied from the 1,875-sequence InterPro-
     confirmed-genuine phaC reference set (PHA_CLEAN_RESULTS.md section
     1.4/4), so the reference structures Foldseek will compare against are
     built from the same verified pool as everything else in this project,
     not an arbitrary external set.
"""
import csv
import pickle
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FA = ROOT / 'PHA_bioprospecting/omdb_search/results'
OUT = Path(__file__).resolve().parent

MIN_LEN, MAX_LEN = 150, 700
N_POSITIVE_CONTROLS = 400
RANDOM_SEED = 42

UNCERTAIN_TIERS = {'No HMM/triad, >=5 other PHA genes', 'No HMM/triad, 1-4 other PHA genes', 'phaC only'}
CONFIDENT_TIERS = {'Catalytic triad complete', 'HMM-supported (no triad)'}

genome_group = pickle.load(open('/tmp/phac_verified_triad_hmm_group.pkl', 'rb'))
genome_targets = pickle.load(open('/tmp/verified_phac_genome_targets.pkl', 'rb'))
target_lengths = pickle.load(open('/tmp/all_phac_target_lengths.pkl', 'rb'))

assignments = {}  # target_id -> cluster representative target_id
with open(FA / 'phaC_cluster0.7_cluster.tsv') as f:
    for row in csv.reader(f, delimiter='\t'):
        assignments[row[1]] = row[0]

# genome metadata for stratifying the positive-control sample by phylum
genome_phylum = {}
with open(FA / 'genome_family_matrix.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        genome_phylum[row['genome']] = row.get('gtdb_phylum') or 'unknown'

def length_ok(cid):
    L = target_lengths.get(cid)
    return L is not None and MIN_LEN <= L <= MAX_LEN

uncertain_clusters = set()
confident_clusters = set()
cluster_phyla = defaultdict(set)  # cluster -> phyla of genomes touching it
for g, tier in genome_group.items():
    phylum = genome_phylum.get(g, 'unknown')
    for t in genome_targets.get(g, []):
        c = assignments.get(t)
        if not c:
            continue
        cluster_phyla[c].add(phylum)
        (uncertain_clusters if tier in UNCERTAIN_TIERS else confident_clusters).add(c)

uncertain_final = sorted(c for c in uncertain_clusters if length_ok(c))
confident_pool = sorted(c for c in confident_clusters if length_ok(c))

# stratified sample of positive controls: proportional to phylum representation
# in the confident pool, so common phyla dominate the sample the same way
# they dominate the real data, rather than an artificial even split
by_phylum = defaultdict(list)
for c in confident_pool:
    phyla = cluster_phyla.get(c, {'unknown'})
    phylum = sorted(phyla)[0]  # deterministic pick if a cluster spans >1 phylum
    by_phylum[phylum].append(c)

random.seed(RANDOM_SEED)
positive_controls = []
remaining = N_POSITIVE_CONTROLS
phyla_sorted = sorted(by_phylum, key=lambda p: -len(by_phylum[p]))
for i, phylum in enumerate(phyla_sorted):
    pool = by_phylum[phylum]
    random.shuffle(pool)
    share = max(1, round(N_POSITIVE_CONTROLS * len(pool) / len(confident_pool)))
    take = pool[:share]
    positive_controls.extend(take)
positive_controls = sorted(set(positive_controls))[:N_POSITIVE_CONTROLS + 50]  # small slack from rounding

print(f'uncertain-tier representatives: {len(uncertain_final)}')
print(f'positive-control representatives: {len(positive_controls)} (from {len(confident_pool)} eligible, {len(by_phylum)} phyla)')

# ---- pull sequences ----
def read_fasta(path):
    cur = None
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur:
                    yield cur, ''.join(seq)
                cur, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip())
        if cur:
            yield cur, ''.join(seq)

wanted = set(uncertain_final) | set(positive_controls)
sequences = {}
for h, s in read_fasta('/tmp/all_phac_dataset.faa'):
    if h in wanted:
        sequences[h] = s

missing = wanted - set(sequences)
if missing:
    print(f'WARNING: {len(missing)} wanted cluster representatives not found in /tmp/all_phac_dataset.faa')

def write_fasta(path, ids):
    n = 0
    with open(path, 'w') as out:
        for cid in ids:
            seq = sequences.get(cid)
            if not seq:
                continue
            out.write(f'>{cid}\n')
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + '\n')
            n += 1
    print(f'wrote {path} ({n} sequences)')

write_fasta(OUT / 'uncertain_cluster_representatives.faa', uncertain_final)
write_fasta(OUT / 'positive_control_representatives.faa', positive_controls)

# ---- reference set: copy from the already-fetched audited FASTA ----
import shutil
ref_src = Path('/tmp/audited_phac_refs.faa')
ref_dst = OUT / 'audited_reference_set.faa'
if ref_src.exists():
    shutil.copy(ref_src, ref_dst)
    n_refs = sum(1 for _ in read_fasta(ref_dst))
    print(f'copied {ref_dst} ({n_refs} sequences) from {ref_src}')
else:
    print(f'WARNING: {ref_src} not found -- re-fetch the audited reference set '
          f'(see figures/PHA_CLEAN_RESULTS.md section 4 for the UniProt bulk-fetch command) before running this again.')

# ---- manifest for the sbatch array job ----
manifest_path = OUT / 'fold_manifest.tsv'
with open(manifest_path, 'w') as out:
    out.write('target_id\tset\n')
    for cid in uncertain_final:
        out.write(f'{cid}\tuncertain\n')
    for cid in positive_controls:
        out.write(f'{cid}\tpositive_control\n')
    for h, _ in read_fasta(ref_dst) if ref_dst.exists() else []:
        out.write(f'{h}\treference\n')
print(f'wrote {manifest_path}')
