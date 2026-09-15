#!/bin/bash
# Combines every batch's hits file (from cluster/run_gopc_search_batched.sbatch's
# job array) into the same <family>_hits.tsv / _unique_targets.tsv /
# _summary.tsv a single-shot gopc-search run would have produced. Fast,
# CPU-only, no GPU needed -- runs on `service`.
#
# Chain it to run automatically once the whole array finishes:
#   sbatch --account=191001-364393 --array=0-3 --export=ALL,FAMILY=phaC,N_BATCHES=4 \
#     cluster/run_gopc_search_batched.sbatch
#   # note the printed job id, e.g. 12345, then:
#   sbatch --account=191001-364393 --dependency=afterok:12345 \
#     --export=ALL,FAMILY=phaC,N_BATCHES=4 cluster/run_gopc_finalize_batches.sbatch
# N_BATCHES here must match what the array job was actually run with, or
# finalize will look for a batch that was never produced and fail loudly.
#SBATCH --job-name=pha_gopc_finalize_batches
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gopc_finalize_batches_%j.out
#SBATCH --error=logs/gopc_finalize_batches_%j.err

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

: "${FAMILY:?Set FAMILY=<family_id>, comma-separated list, or 'all'.}"
: "${N_BATCHES:?Set N_BATCHES=<total batch count> -- must match the array job being finalized.}"
MAX_SEQS="${MAX_SEQS:-10000}"
echo "FAMILY=${FAMILY} N_BATCHES=${N_BATCHES} MAX_SEQS=${MAX_SEQS}"

pha-reference gopc-finalize-batches --family "${FAMILY}" --n-batches "${N_BATCHES}" --max-seqs "${MAX_SEQS}"
pha-reference gopc-combine

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/gopc_search/results/ -- same file layout as an unbatched run."
