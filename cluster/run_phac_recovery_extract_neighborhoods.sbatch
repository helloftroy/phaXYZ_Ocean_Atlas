#!/bin/bash
# Step 2 of 3. Downloads each of the 13,522 qualifying genomes' own gene
# calls directly from OMDB (predictable per-genome URL, confirmed live --
# see phaatlas/pipeline/omdb_metadata.py and phac_recovery.py's docstring),
# locates each genome's known other-pha-family hits by exact sequence
# match against phac_recovery/target_sequences.faa (step 1's output), and
# writes the union of their +/-10-gene neighborhoods to one combined
# FASTA -- so the 13k+ raw downloads don't all need to be kept on disk,
# only the much smaller neighborhood subset that actually gets searched.
#
# NEEDS INTERNET (13,522 individual HTTPS fetches to sunagawalab.ethz.ch)
# -- `service` partition, not GPU, same as
# cluster/run_omdb_enrich_metadata.sbatch. Long walltime because of the
# sheer number of individual requests, not because any one is slow.
#
#   sbatch --account=191001-364393 cluster/run_phac_recovery_extract_neighborhoods.sbatch
# Override the neighborhood window (default +/-10 genes):
#   sbatch --account=191001-364393 --export=ALL,WINDOW=15 cluster/run_phac_recovery_extract_neighborhoods.sbatch
#SBATCH --job-name=phac_recovery_extract_neighborhoods
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --output=logs/phac_recovery_extract_neighborhoods_%j.out
#SBATCH --error=logs/phac_recovery_extract_neighborhoods_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

if [ ! -f cluster/env_activate.sh ]; then
  echo "cluster/env_activate.sh is missing -- run ./cluster/setup_env.sh first." >&2
  exit 2
fi
# shellcheck source=/dev/null
source cluster/env_activate.sh

echo "Job started: $(date)"
echo "Host: $(hostname)"

WINDOW="${WINDOW:-10}"
MANIFEST="phac_recovery/genome_download_manifest.tsv"
HITS="phac_recovery/other_pha_hits.tsv"
TARGET_SEQS="phac_recovery/target_sequences.faa"
OUT_FASTA="phac_recovery/neighborhoods.faa"
OUT_REPORT="phac_recovery/extraction_report.tsv"

for required in "${MANIFEST}" "${HITS}" "${TARGET_SEQS}"; do
  if [ ! -f "${required}" ]; then
    echo "${required} not found -- see cluster/run_phac_recovery_extract_targets.sbatch (step 1) and phac_recovery/build_local_manifests.py." >&2
    exit 2
  fi
done

pha-reference phac-recovery-extract-neighborhoods \
  --genome-download-manifest "${MANIFEST}" \
  --other-pha-hits "${HITS}" \
  --target-sequences-faa "${TARGET_SEQS}" \
  --out-fasta "${OUT_FASTA}" \
  --out-report "${OUT_REPORT}" \
  --window "${WINDOW}"

echo "Job finished: $(date)"
echo "Next: sbatch cluster/run_phac_recovery_hmmsearch.sbatch"
