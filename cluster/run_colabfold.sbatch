#!/bin/bash
# ColabFold structure prediction for the specific candidates ESMFold
# handled poorly (very low pLDDT) or where Foldseek's verdict against the
# reference set was ambiguous -- see PHA_CLEAN_RESULTS.md section 6.3,
# step 4. This is meant for a SMALL escalation subset identified after
# looking at the ESMFold + Foldseek results, not the full ~7,300-sequence
# set (that's what ESMFold was for, precisely because it's cheap enough
# to run at that scale and ColabFold isn't).
#
# Uses the free remote MSA server (api.colabfold.com) by default -- no
# local sequence database needed. Only add --local-msa-dbs (and point it
# at a real local ColabFold/UniRef30 database) if the remote server's
# rate limits become a real bottleneck; not expected to at "small
# escalation subset" scale.
#
# Usage:
#   sbatch --account=191001-364393 \
#     --export=ALL,FASTA=structure_prediction/escalation_candidates.faa,OUT_DIR=structure_prediction/colabfold_out \
#     cluster/run_colabfold.sbatch
#
#SBATCH --job-name=pha_colabfold
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/colabfold_%j.out
#SBATCH --error=logs/colabfold_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs
export PYTHONUNBUFFERED=1

: "${FASTA:?Set FASTA=path to the escalation-subset FASTA (build this from the ESMFold + Foldseek results, not the full input sets).}"
: "${OUT_DIR:?Set OUT_DIR=output directory for ColabFold predictions.}"

SCRATCH_BASE="${SCRATCH_BASE:-/scratch/morrill/users/hmp278}"
CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-${SCRATCH_BASE}/conda_envs/colabfold}"
export COLABFOLD_DATA_DIR="${COLABFOLD_DATA_DIR:-${SCRATCH_BASE}/colabfold_weights}"

MINIFORGE_HOME="${MINIFORGE_HOME:-$HOME/miniforge3}"
# shellcheck source=/dev/null
source "${MINIFORGE_HOME}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PREFIX}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "FASTA=${FASTA} OUT_DIR=${OUT_DIR}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

colabfold_batch "${FASTA}" "${OUT_DIR}"

echo "Job finished: $(date)"
