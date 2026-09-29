"""Deeper follow-up on the Scalindua/Pseudothioglobus identical-phaC finding
(§12.6): geographic distance between the two genomes, rarity of their
shared 100%-identity variant within its phaC_cluster0.7 group (compared
against a second cross-genus pair in the same cluster, as a contrast
case), and the full per-family gene detail for the one Scalindua genome
that carries phaC. All from primary tables -- no network fetch needed
(completeness/contamination/lat-lon are already cached in
phaC_unique_targets_with_metadata.tsv / genome_family_matrix.tsv).
"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
FA = ROOT / 'PHA_bioprospecting' / 'omdb_search' / 'results'

REP = 'OMDBv2.0_AA_G_NR100_000013537294'
SCALINDUA_GENOME = 'BEMA21-1_SAMN15000289_MAG_00000040'
PSEUDOTHIOGLOBUS_GENOME = 'GLAS15-1_SAMN02905564_MAG_00000008'
FAMILIES = ['phaA', 'phaB', 'phaC', 'phaD', 'phaE', 'phaF', 'phaG', 'phaI',
            'phaJ', 'phaP', 'phaQ', 'phaR_regulator', 'phaR_synthase', 'phaY', 'phaZ']


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


# 1. cluster0.7 membership for REP's cluster
member_to_rep = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) < 2:
            continue
        if row[0] == REP:
            member_to_rep[row[1]] = row[0]
members = set(member_to_rep)

rows = []
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in members:
            rows.append(row)

from collections import defaultdict, Counter
by_target = defaultdict(list)
for r in rows:
    by_target[r['target_id']].append(r)

variants_out = ROOT / 'figures' / 'scalindua_pseudothioglobus_cluster_variants.tsv'
with open(variants_out, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['target_id', 'n_genomes', 'genera', 'is_our_pair'])
    for t, rs in sorted(by_target.items(), key=lambda kv: -len(kv[1])):
        genera = ','.join(f'{g}:{n}' for g, n in Counter(r['gtdb_genus'] for r in rs).most_common())
        writer.writerow([t, len(rs), genera, t == REP])
print(f'Wrote {variants_out} ({len(by_target)} distinct 100%-identity variants in this cluster)')

# 2. geographic distance for our pair, and the comparison cross-genus pair
def genome_latlon(genome):
    with open(FA / 'genome_family_matrix.tsv', newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['genome'] == genome:
                return float(row['latitude_degN']), float(row['longitude_degE']), row['study_id']
    return None

lat1, lon1, study1 = genome_latlon(SCALINDUA_GENOME)
lat2, lon2, study2 = genome_latlon(PSEUDOTHIOGLOBUS_GENOME)
d = haversine_km(lat1, lon1, lat2, lon2)
print(f'{SCALINDUA_GENOME} ({study1}) <-> {PSEUDOTHIOGLOBUS_GENOME} ({study2}): {d:.0f} km apart')

# find the OTHER cross-genus pair in this cluster (contrast case)
for t, rs in by_target.items():
    genera = set(r['gtdb_genus'] for r in rs)
    if len(rs) == 2 and len(genera) == 2 and t != REP:
        (g1, g2) = rs
        d2 = haversine_km(float(g1['latitude_degN']), float(g1['longitude_degE']),
                           float(g2['latitude_degN']), float(g2['longitude_degE']))
        print(f'Contrast pair {t}: {g1["genome"]} ({g1["study_id"]}) <-> {g2["genome"]} ({g2["study_id"]}): '
              f'{d2:.0f} km apart, same study={g1["study_id"] == g2["study_id"]}')

# 3. full per-family gene detail for the one Scalindua phaC-positive genome
genes_out = ROOT / 'figures' / 'scalindua_bema21_other_pha_genes.tsv'
gene_rows = []
for fam in FAMILIES:
    path = FA / f'{fam}_unique_targets_with_metadata.tsv'
    if not path.exists():
        continue
    with open(path, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('genome') == SCALINDUA_GENOME:
                gene_rows.append({
                    'family': fam, 'target_id': row['target_id'], 'best_query': row['best_query'],
                    'best_pident': row['best_pident'], 'best_qcov': row['best_qcov'], 'best_tcov': row['best_tcov'],
                })
with open(genes_out, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=['family', 'target_id', 'best_query', 'best_pident', 'best_qcov', 'best_tcov'])
    writer.writeheader()
    writer.writerows(gene_rows)
print(f'Wrote {genes_out} ({len(gene_rows)} family hits for {SCALINDUA_GENOME})')
