#!/bin/bash
# Pulls out the actual protein sequences for a family's discovered
# target_ids from the mmseqs target database ALREADY BUILT for the search
# stage (run_omdb_build_db_gpu.sbatch) -- no need to re-download
# OMDBv2.0_AA_G_NR100.faa.gz (~46GB). Uses `mmseqs createsubdb --id-mode 1`
# (select by FASTA identifier rather than internal DB key -- confirmed
# live against a small local test target db) then `mmseqs convert2fasta`.
# Neither needs a GPU, so this runs on `service` with the plain CPU
# mmseqs build (cluster/bin/mmseqs) even though the target db itself was
# built GPU-compatible -- that storage format is still a normal readable
# sequence DB for non-search operations.
#
# This is step 1 of "cluster PhaC more deeply than the best_query proxy":
#   1. THIS script -- extract sequences for a family's target_ids
#   2. cluster/run_omdb_cluster_sequences.sbatch -- cluster them at a real identity threshold
#   3. `pha-reference cluster-ecology` / `sequence-embedding` (local, fast, no scheduling needed)
#
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_extract_cluster_sequences.sbatch
#SBATCH --job-name=pha_omdb_extract_cluster_seqs
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_extract_cluster_seqs_%j.out
#SBATCH --error=logs/omdb_extract_cluster_seqs_%j.err

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

: "${FAMILY:?Set FAMILY to a single family_id.}"
TARGET_DB="${TARGET_DB:-PHA_bioprospecting/omdb_search/target_db_gpu/omdb_db}"
MMSEQS_BIN="${MMSEQS_BIN:-$(pwd)/cluster/bin/mmseqs}"
RESULTS_DIR="${RESULTS_DIR:-PHA_bioprospecting/omdb_search/results}"
echo "FAMILY=${FAMILY} TARGET_DB=${TARGET_DB}"

if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build it first (see cluster/run_omdb_build_db_gpu.sbatch)." >&2
  exit 2
fi
if [ ! -x "${MMSEQS_BIN}" ]; then
  echo "${MMSEQS_BIN} not found/executable -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi

pha-reference cluster-extract-ids --family "${FAMILY}" --results-dir "${RESULTS_DIR}"

WANTED_IDS="${RESULTS_DIR}/${FAMILY}_cluster_wanted_ids.txt"
SUBSET_DB="${RESULTS_DIR}/${FAMILY}_cluster_seq_db/${FAMILY}_subset"
FASTA_OUT="${RESULTS_DIR}/${FAMILY}_cluster_sequences.faa"
mkdir -p "$(dirname "${SUBSET_DB}")"

"${MMSEQS_BIN}" createsubdb --id-mode 1 "${WANTED_IDS}" "${TARGET_DB}" "${SUBSET_DB}"
"${MMSEQS_BIN}" convert2fasta "${SUBSET_DB}" "${FASTA_OUT}"

echo "Job finished: $(date)"
echo "Extracted sequences: ${FASTA_OUT}"
grep -c "^>" "${FASTA_OUT}" | xargs echo "Sequence count:"
