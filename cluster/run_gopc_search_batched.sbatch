#!/bin/bash
# Batched GPU search: for a family whose full query set doesn't finish
# against GOPC inside one 48h job, splits its queries into N_BATCHES
# independent mmseqs runs (a SLURM job array, one task per batch) so each
# individual run is short enough to actually complete, rather than
# resubmitting the full-family search over and over.
#
# WHY THIS EXISTS (confirmed live): a killed/timed-out
# cluster/run_gopc_search.sbatch run CANNOT be resumed by resubmitting the
# same command -- mmseqs creates a fresh, randomly-named tmp subdirectory
# on every invocation, so a resubmit restarts the entire ~2.46-billion-
# sequence GOPC scan from zero every time. Splitting the QUERY side into
# batches does NOT reduce the total GPU-time needed (each batch still
# scans all of GOPC once, so total cost goes UP relative to one
# hypothetical uninterrupted run) -- what it buys is turning "one run that
# might never finish inside the time limit" into several independent runs
# that each definitely will, at the cost of more total GPU-hours.
#
# --array MUST cover exactly 0..N_BATCHES-1 -- e.g. for N_BATCHES=4:
#   sbatch --account=191001-364393 --array=0-3 --export=ALL,FAMILY=phaC,N_BATCHES=4 \
#     cluster/run_gopc_search_batched.sbatch
# Once every array task finishes, finalize (combines the batch hits files
# into the same <family>_hits.tsv/_unique_targets.tsv/_summary.tsv a
# single-shot run would produce) -- chain it with --dependency so it
# starts automatically once the whole array completes:
#   sbatch --account=191001-364393 --dependency=afterok:<array_job_id> \
#     --export=ALL,FAMILY=phaC,N_BATCHES=4 cluster/run_gopc_finalize_batches.sbatch
#
# How many batches: no formula here (we don't have solid enough timing
# data yet to derive one) -- pick N_BATCHES so that (this family's query
# count / N_BATCHES) feels comparable to a family that's already finished
# successfully in one 48h shot, and adjust from what you observe. Each
# batch's own log shows "Query database size: N" near the top, same as an
# unbatched run.
#
# Same --mem/--cpus-per-task/--gres reasoning as run_gopc_search.sbatch --
# each batch task still scans the full target, so needs the same resources
# as a full run would, just for fewer queries.
#SBATCH --job-name=pha_gopc_search_batch
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gopc_search_batch_%A_%a.out
#SBATCH --error=logs/gopc_search_batch_%A_%a.err

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

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/gopc_search/target_db_gpu/gopc_db}"
if [ ! -f "${TARGET_DB}" ]; then
  echo "${TARGET_DB} not found -- build it first: cluster/run_gopc_build_db_gpu.sbatch" >&2
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
pha-reference gopc-search-batch \
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
echo "Once every batch has completed, run cluster/run_gopc_finalize_batches.sbatch."
