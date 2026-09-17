"""Companion export for figures/phaC_regional_specialists.png -- one CSV
row per selected cluster: identifying/ecological info + the cluster's own
representative sequence. The mmseqs cluster_id IS its representative
member's own sequence ID (see phylogenetics.py's module docstring), so
this is that member's real sequence, not a re-aligned/re-derived consensus.

Cluster list, per-cluster genus/depth/position values, and the "drop
missing depths before averaging" convention are kept identical to
plot_regional_specialists.py so this CSV and that figure always describe
the same 12 clusters -- if that list changes, update both files.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from phaatlas.pipeline import sequence_clustering as sc

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = Path(__file__).resolve().parent.parent

SELECTED = [
    '000033545441', '000034200129', '000140439260', '000067604172',
    '000045952813', '000005303736', '000232749406', '000012716247',
    '000226102035', '000018400463', '000001867111', '000027486053',
]


def read_fasta(path):
    cur, seq = None, []
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


eco_rows = {r['cluster_id'][-12:]: r for r in csv.DictReader(open(FA / 'phaC_cluster0.7_cluster_ecology.tsv', newline=''), delimiter='\t')}
full_ids = {short: eco_rows[short]['cluster_id'] for short in SELECTED}

sequences = {}
for hid, seq in read_fasta(FA / 'phaC_all_representatives.faa'):
    if hid in full_ids.values():
        sequences[hid] = seq

assignments = sc.load_cluster_assignments(FA / 'phaC_cluster0.7_cluster.tsv')
depths = {fid: [] for fid in full_ids.values()}
seen_genomes = {fid: set() for fid in full_ids.values()}
with open(FA / 'phaC_unique_targets_with_metadata_depth.tsv', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        cid = assignments.get(row.get('target_id', ''))
        if cid not in depths:
            continue
        genome = row.get('genome', '')
        if not genome or genome in seen_genomes[cid]:
            continue
        seen_genomes[cid].add(genome)
        d = row.get('depth_m', '')
        if d:
            depths[cid].append(float(d))

out_path = OUT / 'phaC_regional_specialists_sequences.csv'
with open(out_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['cluster_id', 'dominant_genus', 'n_genomes', 'n_genomes_with_depth', 'mean_depth_m',
                      'latitude_degN', 'longitude_degE', 'sequence'])
    for short in SELECTED:
        fid = full_ids[short]
        r = eco_rows[short]
        ds = depths[fid]
        mean_depth = round(sum(ds) / len(ds), 2) if ds else ''
        writer.writerow([
            fid, r['top_genera'].split(' (')[0], r['n_genomes'], len(ds), mean_depth,
            r['geo_centroid_lat'], r['geo_centroid_lon'], sequences[fid],
        ])

print('saved', out_path)
print(f'{len(sequences)} sequences written')
