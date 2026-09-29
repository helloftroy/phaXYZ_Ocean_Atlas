"""Per-genome phaC breakdown for Methylocystis/Methylosinus (Type II
methanotrophs, Alphaproteobacteria/Beijerinckiaceae) -- built independently
from this project's own primary tables (genome_family_matrix.tsv,
phaC_unique_targets_with_metadata.tsv, the QC exclusion list, the
alignment-column triad table, and phaC_cluster0.7_cluster.tsv), not from
any other agent's pre-aggregated summary files. Mirrors the structure of
the existing Thioglobus/SUP05 per-genome breakdown
(figures/thioglobus_sup05_phac_taxonomic_breakdown.tsv) so the two are
directly comparable, but is a fresh build: this project's own convention
(see PHA_CLEAN_RESULTS.md) is to independently re-derive rather than trust
another agent's aggregate numbers wholesale, since at least one such prior
result (a Scalindua architecture read) was found to have missed real
pha-family gene calls.

Cross-check: this script's own 34 QC-passing target rows / 20 phaC-positive
genomes / 26 triad-complete rows / 19 genomes-with-triad-complete
reproduces the same four numbers already on record for this group --
independent confirmation that this specific tally is solid, whatever else
in that other pass needs re-checking.
"""
import csv
import sys
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / 'figures' / 'scripts'))
sys.path.insert(0, str(ROOT / 'catalytic_domain'))
import _phac_qc
from find_structural_triad import structural_triad_complete_batch

FA = ROOT / 'PHA_bioprospecting' / 'omdb_search' / 'results'
GENUS_OK = {'Methylocystis', 'Methylosinus'}
OUT = ROOT / 'figures' / 'methylocystis_phac_taxonomic_breakdown.tsv'

# 1. phaC-positive genomes of these two genera, genome-level call
phac_pos_genomes = set()
species_by_genome = {}
with open(FA / 'genome_family_matrix.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['gtdb_genus'] in GENUS_OK and int(row['n_phaC']) > 0:
            phac_pos_genomes.add(row['genome'])
            species_by_genome[row['genome']] = (row['gtdb_genus'], row['gtdb_species'])
print(f'{len(phac_pos_genomes)} phaC-positive Methylocystis/Methylosinus genomes')

# 2. raw phaC target rows for these genomes, QC-filtered
raw_rows = []
with open(FA / 'phaC_unique_targets_with_metadata.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['genome'] in phac_pos_genomes:
            raw_rows.append(row)
bad = _phac_qc.load_bad_targets()
qc_rows = [r for r in raw_rows if r['target_id'] not in bad]
print(f'{len(raw_rows)} raw target rows -> {len(qc_rows)} QC-passing '
      f'({len(set(r["genome"] for r in qc_rows))} genomes)')

# 3. alignment-column triad + structural triad (where a local PDB exists)
target_ids = sorted(set(r['target_id'] for r in qc_rows))
triad_alignment = {}
with open(ROOT / 'catalytic_domain' / 'phac_catalytic_triad.tsv', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['target_id'] in set(target_ids):
            triad_alignment[row['target_id']] = row['triad_complete'] == 'True'
triad_structural = structural_triad_complete_batch(target_ids)
print(f'{sum(triad_alignment.values())}/{len(target_ids)} distinct targets triad-complete (alignment-column); '
      f'{len(triad_structural)} have a local PDB, {sum(triad_structural.values())} of those triad-complete structurally')

n_triad_true_rows = sum(1 for r in qc_rows if triad_alignment.get(r['target_id'], False))
genomes_with_triad = set(r['genome'] for r in qc_rows if triad_alignment.get(r['target_id'], False))
print(f'{n_triad_true_rows} QC-passing target-assignment rows triad-complete; '
      f'{len(genomes_with_triad)} genomes with >=1 triad-complete target')

# 4. cluster0.7 assignment
member_to_rep = {}
with open(FA / 'phaC_cluster0.7_cluster.tsv', newline='') as f:
    for row in csv.reader(f, delimiter='\t'):
        if len(row) < 2:
            continue
        member_to_rep[row[1]] = row[0]

# 5. write per-genome breakdown
by_genome = defaultdict(list)
for r in qc_rows:
    by_genome[r['genome']].append(r)

fields = [
    'genome', 'gtdb_genus', 'gtdb_species', 'n_qc_targets', 'target_ids',
    'cluster0.7_representatives', 'any_triad_complete_alignment',
    'any_triad_complete_structural', 'has_local_pdb', 'completeness', 'contamination',
]
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
    writer.writeheader()
    for g in sorted(by_genome):
        rows = by_genome[g]
        genus, species = species_by_genome.get(g, ('', ''))
        tids = sorted(r['target_id'] for r in rows)
        reps = sorted({member_to_rep.get(t, '') for t in tids})
        writer.writerow({
            'genome': g,
            'gtdb_genus': genus,
            'gtdb_species': species,
            'n_qc_targets': len(rows),
            'target_ids': ','.join(tids),
            'cluster0.7_representatives': ','.join(reps),
            'any_triad_complete_alignment': any(triad_alignment.get(t, False) for t in tids),
            'any_triad_complete_structural': any(triad_structural.get(t, False) for t in tids),
            'has_local_pdb': any(t in triad_structural for t in tids),
            'completeness': rows[0]['completeness'],
            'contamination': rows[0]['contamination'],
        })

print(f'Wrote {len(by_genome)} genomes -> {OUT}')
