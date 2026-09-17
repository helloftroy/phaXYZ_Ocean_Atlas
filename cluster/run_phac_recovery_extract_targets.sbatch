#!/bin/bash
# Step 1 of 3 for "find the missed phaC in genomes with >5 other PHA genes
# but no phaC hit" (see phaatlas/pipeline/phac_recovery.py's module
# docstring for the full rationale). Pulls the actual sequences for the
# ~170k distinct target_ids that are the OTHER-pha-family anchors in our
# 13,522 qualifying genomes, from the mmseqs target_db ALREADY BUILT for
# the original search (same trick as
# cluster/run_omdb_extract_cluster_sequences.sbatch: `mmseqs createsubdb
# --id-mode 1` + `convert2fasta`, no need to touch the 49GB
# OMDBv2.0_AA_G_NR100.faa.gz or the 4.2GB cluster.tsv at all).
#
# Needs phac_recovery/target_ids.txt on this machine first -- scp it up
# from the fair_ocean_agent/phac_recovery/ output of
# phac_recovery/build_local_manifests.py (run locally, not on the
# cluster -- it only needs files already sitting in fair_ocean_agent/):
#   scp fair_ocean_agent/phac_recovery/target_ids.txt <cluster>:<repo>/phac_recovery/
#   scp fair_ocean_agent/phac_recovery/other_pha_hits.tsv <cluster>:<repo>/phac_recovery/
#   scp fair_ocean_agent/phac_recovery/genome_download_manifest.tsv <cluster>:<repo>/phac_recovery/
#
#   sbatch --account=191001-364393 cluster/run_phac_recovery_extract_targets.sbatch
#SBATCH --job-name=phac_recovery_extract_targets
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/phac_recovery_extract_targets_%j.out
#SBATCH --error=logs/phac_recovery_extract_targets_%j.err

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

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/omdb_search/target_db_gpu/omdb_db}"
MMSEQS_BIN="${MMSEQS_BIN:-$(pwd)/cluster/bin/mmseqs}"
TARGET_IDS="phac_recovery/target_ids.txt"
SUBSET_DB="phac_recovery/target_db_subset/phac_recovery_subset"
FASTA_OUT="phac_recovery/target_sequences.faa"

if [ ! -f "${TARGET_IDS}" ]; then
  echo "${TARGET_IDS} not found -- scp it from fair_ocean_agent/phac_recovery/ (built by phac_recovery/build_local_manifests.py locally)." >&2
  exit 2
fi
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build it first (see cluster/run_omdb_build_db_gpu.sbatch)." >&2
  exit 2
fi
if [ ! -x "${MMSEQS_BIN}" ]; then
  echo "${MMSEQS_BIN} not found/executable -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi

mkdir -p "$(dirname "${SUBSET_DB}")"
N_IDS=$(wc -l < "${TARGET_IDS}")
echo "Extracting sequences for ${N_IDS} target_ids..."

"${MMSEQS_BIN}" createsubdb --id-mode 1 "${TARGET_IDS}" "${TARGET_DB}" "${SUBSET_DB}"
"${MMSEQS_BIN}" convert2fasta "${SUBSET_DB}" "${FASTA_OUT}"

echo "Job finished: $(date)"
echo "Extracted sequences: ${FASTA_OUT}"
grep -c "^>" "${FASTA_OUT}" | xargs echo "Sequence count:"
echo "Next: sbatch cluster/run_phac_recovery_extract_neighborhoods.sbatch"
