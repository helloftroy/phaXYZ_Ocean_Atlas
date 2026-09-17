"""Phase 1: everything needed to hand off to the cluster for the actual
phaC-recovery search, built from genome_family_matrix.tsv and the 14
<family>_unique_targets_with_metadata.tsv files (wherever those actually
are on THIS machine -- see --data-dir below) plus the small OMDB catalog
(fetched fresh here rather than committed, since it's a public,
re-fetchable download, not a search result).

Where the input data lives differs by machine: on the author's Mac it's
fair_ocean_agent/ (a local scp landing zone, not part of this repo); on
the cluster it's more likely PHA_bioprospecting/omdb_search/results/
(wherever `pha-reference pathway-architecture` was actually run against).
Pass --data-dir explicitly if the default guess below is wrong for your
setup -- the error messages below will say clearly if a required file
isn't found there.

Outputs always go to phac_recovery/ next to this script (i.e. the repo
root's phac_recovery/, matching where the sbatch scripts and the
already-committed HMM files expect them), regardless of --data-dir:
  - qualifying_genomes.txt            one genome id per line
  - other_pha_hits.tsv                genome, family, target_id (anchors)
  - genome_download_manifest.tsv      genome, genes_aa_url
  - target_ids.txt                    distinct target_ids needing sequence
                                       extraction from the cluster's
                                       target_db (see
                                       cluster/run_phac_recovery_extract_targets.sbatch)

Usage:
  .venv/bin/python phac_recovery/build_local_manifests.py --data-dir /path/to/results
"""
import argparse
import csv
import sys
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from phaatlas.pipeline import phac_recovery as pr

DEFAULT_DATA_DIR = REPO_ROOT / 'PHA_bioprospecting' / 'omdb_search' / 'results'

parser = argparse.ArgumentParser()
parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA_DIR,
                     help='directory containing genome_family_matrix.tsv and the '
                          '<family>_unique_targets_with_metadata.tsv files '
                          f'(default: {DEFAULT_DATA_DIR})')
args = parser.parse_args()

DATA_DIR = args.data_dir
OUT = Path(__file__).resolve().parent  # phac_recovery/ at repo root, NOT data-dir-relative
OUT.mkdir(parents=True, exist_ok=True)

matrix_path = DATA_DIR / 'genome_family_matrix.tsv'
if not matrix_path.exists():
    sys.exit(f'{matrix_path} not found. Pass --data-dir pointing at wherever '
              f'genome_family_matrix.tsv actually is (run `pha-reference pathway-architecture` '
              f'first if it does not exist yet anywhere).')

CATALOG_URL = 'https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs/OMDBv2.0_data.tsv.gz'
CATALOG_PATH = OUT / 'OMDBv2.0_data.tsv'

FAMILIES = ['phaA', 'phaB', 'phaD', 'phaE', 'phaF', 'phaG', 'phaI', 'phaJ',
            'phaP', 'phaQ', 'phaR_regulator', 'phaR_synthase', 'phaY', 'phaZ']
# phaM excluded: confirmed earlier (session data) that 0 genomes anywhere
# have a phaM hit, so phaM_unique_targets_with_metadata.tsv doesn't exist
# to scan.

print(f'Reading input data from {DATA_DIR}')
qualifying = pr.identify_qualifying_genomes(matrix_path, min_other_genes=5)
print(f'{len(qualifying)} qualifying genomes (no phaC, >5 other PHA genes)')

with open(OUT / 'qualifying_genomes.txt', 'w') as f:
    for g in sorted(qualifying):
        f.write(g + '\n')

family_paths = {fam: DATA_DIR / f'{fam}_unique_targets_with_metadata.tsv' for fam in FAMILIES}
missing = [fam for fam, p in family_paths.items() if not p.exists()]
if missing:
    print(f'WARNING: missing metadata files for {missing} in {DATA_DIR}, skipping those families')
    family_paths = {fam: p for fam, p in family_paths.items() if fam not in missing}
if not family_paths:
    sys.exit(f'None of the 14 <family>_unique_targets_with_metadata.tsv files were found in {DATA_DIR} -- nothing to anchor on.')

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
