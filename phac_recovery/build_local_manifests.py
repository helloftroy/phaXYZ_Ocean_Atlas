"""Phase 1 (local): everything needed to hand off to the cluster for the
actual phaC-recovery search, built from data already sitting in
fair_ocean_agent/ plus the small OMDB catalog (fetched fresh here rather
than committed, since it's a public, re-fetchable download, not a search
result).

Outputs (written to fair_ocean_agent/phac_recovery/, the usual scp
landing zone for cluster-scale artifacts, not the git repo):
  - qualifying_genomes.txt            one genome id per line
  - other_pha_hits.tsv                genome, family, target_id (anchors)
  - genome_download_manifest.tsv      genome, genes_aa_url
  - target_ids.txt                    distinct target_ids needing sequence
                                       extraction from the cluster's
                                       target_db (see
                                       cluster/run_phac_recovery_extract_targets.sbatch)
"""
import csv
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from phaatlas.pipeline import phac_recovery as pr

FA = Path('/Users/hellpark/multimodal_seusmbol/fair_ocean_agent')
OUT = FA / 'phac_recovery'
OUT.mkdir(parents=True, exist_ok=True)

CATALOG_URL = 'https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs/OMDBv2.0_data.tsv.gz'
CATALOG_PATH = OUT / 'OMDBv2.0_data.tsv'

FAMILIES = ['phaA', 'phaB', 'phaD', 'phaE', 'phaF', 'phaG', 'phaI', 'phaJ',
            'phaP', 'phaQ', 'phaR_regulator', 'phaR_synthase', 'phaY', 'phaZ']
# phaM excluded: confirmed earlier (session data) that 0 genomes anywhere
# have a phaM hit, so phaM_unique_targets_with_metadata.tsv doesn't exist
# to scan.

qualifying = pr.identify_qualifying_genomes(FA / 'genome_family_matrix.tsv', min_other_genes=5)
print(f'{len(qualifying)} qualifying genomes (no phaC, >5 other PHA genes)')

with open(OUT / 'qualifying_genomes.txt', 'w') as f:
    for g in sorted(qualifying):
        f.write(g + '\n')

family_paths = {fam: FA / f'{fam}_unique_targets_with_metadata.tsv' for fam in FAMILIES}
missing = [fam for fam, p in family_paths.items() if not p.exists()]
if missing:
    print(f'WARNING: missing metadata files for {missing}, skipping those families')
    family_paths = {fam: p for fam, p in family_paths.items() if fam not in missing}

hits = pr.collect_other_pha_hits(qualifying, family_paths)
print(f'{len(hits)} (genome, family, target_id) anchor rows collected')

with open(OUT / 'other_pha_hits.tsv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['genome', 'family', 'target_id'])
    for h in hits:
        writer.writerow([h.genome, h.family, h.target_id])

distinct_target_ids = sorted({h.target_id for h in hits})
with open(OUT / 'target_ids.txt', 'w') as f:
    for t in distinct_target_ids:
        f.write(t + '\n')
print(f'{len(distinct_target_ids)} distinct target_ids need sequence extraction from the cluster target_db')

if not CATALOG_PATH.exists():
    print(f'Fetching OMDB catalog ({CATALOG_URL})...')
    import gzip
    gz_path = OUT / 'OMDBv2.0_data.tsv.gz'
    urllib.request.urlretrieve(CATALOG_URL, gz_path)
    with gzip.open(gz_path, 'rt') as f_in, open(CATALOG_PATH, 'w') as f_out:
        f_out.write(f_in.read())

n = pr.write_genome_download_manifest(qualifying, CATALOG_PATH, OUT / 'genome_download_manifest.tsv')
print(f'{n} genomes matched in OMDB catalog -> genome_download_manifest.tsv')
if n != len(qualifying):
    print(f'WARNING: {len(qualifying) - n} qualifying genomes NOT found in the OMDB catalog')
