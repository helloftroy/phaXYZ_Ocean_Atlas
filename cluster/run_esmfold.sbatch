#!/bin/bash
# ESMFold structure prediction, batched array job. Run once per input
# FASTA (uncertain_cluster_representatives.faa, positive_control_
# representatives.faa, audited_reference_set.faa -- see
# structure_prediction/build_fold_inputs.py), pointing FASTA/OUT_DIR at
# each in turn.
#
# --array MUST cover exactly 0..N_BATCHES-1. Sizing guide (very rough,
# no live timing confirmed yet -- see structure_prediction/run_esmfold.py's
# own note that this hasn't been run end-to-end on a real GPU node):
#   uncertain_cluster_representatives.faa  (4,992 seqs, median 413aa) -> try N_BATCHES=40
#   positive_control_representatives.faa     (428 seqs, median 420aa) -> try N_BATCHES=4
#   audited_reference_set.faa              (1,875 seqs, median 578aa) -> try N_BATCHES=20
# Adjust after watching how long the first few array tasks actually take
# (each task's own stdout logs a per-sequence time) -- these are starting
# guesses sized to keep each task comfortably under a few hours, not a
# measured throughput number.
#
# Usage:
#   sbatch --account=191001-364393 --array=0-39 \
#     --export=ALL,FASTA=structure_prediction/uncertain_cluster_representatives.faa,OUT_DIR=structure_prediction/esmfold_out/uncertain,N_BATCHES=40 \
#     cluster/run_esmfold.sbatch
#
# Resumable: run_esmfold.py skips any target_id whose .pdb already exists
# in OUT_DIR, so re-submitting the same array after a partial failure (OOM
# on a few long sequences, a node preemption, etc.) only redoes the gaps.
#
#SBATCH --job-name=pha_esmfold
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/esmfold_%A_%a.out
#SBATCH --error=logs/esmfold_%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset -- submit with --array=0-<N_BATCHES-1>, not a plain sbatch.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

: "${FASTA:?Set FASTA=<path to one of the three structure_prediction/*.faa input files>.}"
: "${OUT_DIR:?Set OUT_DIR=<output directory for this FASTA's PDB files>.}"
: "${N_BATCHES:?Set N_BATCHES=<total batch count> -- must match the --array range used (0..N_BATCHES-1).}"

MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-esmfold}"
# shellcheck source=/dev/null
source "${MINIFORGE_HOME}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Array task: ${SLURM_ARRAY_TASK_ID}/${N_BATCHES}"
echo "FASTA=${FASTA} OUT_DIR=${OUT_DIR}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

python structure_prediction/run_esmfold.py \
  --fasta "${FASTA}" \
  --out-dir "${OUT_DIR}" \
  --n-batches "${N_BATCHES}" \
  --batch-index "${SLURM_ARRAY_TASK_ID}"

echo "Job finished: $(date)"
echo "Batch ${SLURM_ARRAY_TASK_ID}/${N_BATCHES} done for ${FASTA}."
