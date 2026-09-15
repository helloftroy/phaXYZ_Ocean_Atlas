#!/bin/bash
# Batched GPU search against OMDB, for parity with
# cluster/run_gopc_search_batched.sbatch (see that script for the full
# rationale) -- included in case a family turns out to need it against
# OMDB too, even though OMDB's ~10x smaller size makes that less likely
# than it was for GOPC.
#
# --array MUST cover exactly 0..N_BATCHES-1 -- e.g. for N_BATCHES=4:
#   sbatch --account=191001-364393 --array=0-3 --export=ALL,FAMILY=phaC,N_BATCHES=4 \
#     cluster/run_omdb_search_batched.sbatch
# Then finalize once the whole array completes:
#   sbatch --account=191001-364393 --dependency=afterok:<array_job_id> \
#     --export=ALL,FAMILY=phaC,N_BATCHES=4 cluster/run_omdb_finalize_batches.sbatch
#
# --mem=120G: matches run_gopc_search.sbatch's own value -- confirmed live
# OOM-killed in convertalis (final output-writing step) at 100G for a GOPC
# family with a large result set, see that script's header for detail.
#SBATCH --job-name=pha_omdb_search_batch
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_search_batch_%A_%a.out
#SBATCH --error=logs/omdb_search_batch_%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset -- submit with --array=0-<N_BATCHES-1>, not a plain sbatch.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

if [ ! -f cluster/env_activate.sh ]; then
  echo "cluster/env_activate.sh is missing -- run ./cluster/setup_env.sh first." >&2
  exit 2
fi
# shellcheck source=/dev/null
source cluster/env_activate.sh

MMSEQS_GPU_BIN="${MMSEQS_GPU_BIN:-$(pwd)/cluster/bin-gpu/mmseqs}"
if [ ! -x "${MMSEQS_GPU_BIN}" ]; then
  echo "${MMSEQS_GPU_BIN} not found/executable -- install it first:" >&2
  echo "  MMSEQS_VARIANT=linux-gpu MMSEQS_OUTPUT_DIR=cluster/bin-gpu ./cluster/install_mmseqs2.sh" >&2
  exit 2
fi
echo "mmseqs (GPU build): $("${MMSEQS_GPU_BIN}" version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/omdb_search/target_db_gpu/omdb_db}"
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build it first: cluster/run_omdb_build_db_gpu.sbatch" >&2
  exit 2
fi

: "${FAMILY:?Set FAMILY=<family_id> (single family, not comma-separated/'all').}"
: "${N_BATCHES:?Set N_BATCHES=<total batch count> -- must match the --array range used (0..N_BATCHES-1).}"
SENSITIVITY="${SENSITIVITY:-7.0}"
EVALUE="${EVALUE:-1e-5}"
COVERAGE="${COVERAGE:-0.0}"
MAX_SEQS="${MAX_SEQS:-10000}"
echo "FAMILY=${FAMILY} N_BATCHES=${N_BATCHES} BATCH_INDEX=${SLURM_ARRAY_TASK_ID} MAX_SEQS=${MAX_SEQS}"

pha-reference init-db
pha-reference export-queries
pha-reference omdb-search-batch \
  --family "${FAMILY}" \
  --n-batches "${N_BATCHES}" \
  --batch-index "${SLURM_ARRAY_TASK_ID}" \
  --target-db "${TARGET_DB}" \
  --sensitivity "${SENSITIVITY}" \
  --evalue "${EVALUE}" \
  --coverage "${COVERAGE}" \
  --max-seqs "${MAX_SEQS}" \
  --threads "${SLURM_CPUS_PER_TASK:-16}" \
  --mmseqs-bin "${MMSEQS_GPU_BIN}" \
  --gpu

echo "Job finished: $(date)"
echo "Batch ${SLURM_ARRAY_TASK_ID}/${N_BATCHES} done for ${FAMILY}."
echo "Once every batch has completed, run cluster/run_omdb_finalize_batches.sbatch."
