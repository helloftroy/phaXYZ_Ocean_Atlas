#!/bin/bash
# Structure-comparison step from PHA_CLEAN_RESULTS.md section 6.3/9.7:
# Foldseek search of query structures (uncertain/positive_control/the new
# all_phac_dedup set) against the folded, audited reference structures.
# This is the sbatch version of structure_prediction/run_foldseek_search.sh
# -- same logic, moved off an interactive terminal since a full run took
# noticeably longer than the ESMFold step and is not safe to leave
# attached to an open SSH session (confirmed live: it ran long enough to
# risk being killed by a dropped connection).
#
# Single job, not an array: unlike ESMFold (one GPU, one sequence at a
# time, so more array tasks = more parallel GPUs), Foldseek's own search
# is already multi-threaded across --cpus-per-task for the whole query
# set at once -- splitting it into array tasks would not help and would
# just rebuild the reference database redundantly in each task.
#
# Usage (after ESMFold has produced at least the reference + one query
# set's .pdb files under structure_prediction/esmfold_out/):
#   sbatch --account=191001-364393 \
#     --export=ALL,QUERY_SETS="uncertain positive_control" \
#     cluster/run_foldseek.sbatch
#   # or, for the comprehensive all_phac_dedup run once it has been folded:
#   sbatch --account=191001-364393 \
#     --export=ALL,QUERY_SETS=all_phac_dedup \
#     cluster/run_foldseek.sbatch
#
#SBATCH --job-name=pha_foldseek
#SBATCH --partition=service
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/foldseek_%j.out
#SBATCH --error=logs/foldseek_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs

export FOLDSEEK_BIN="${FOLDSEEK_BIN:-cluster/bin/foldseek}"
export ESMFOLD_OUT="${ESMFOLD_OUT:-structure_prediction/esmfold_out}"
export OUT_DIR="${OUT_DIR:-structure_prediction/foldseek_out}"
export QUERY_SETS="${QUERY_SETS:-uncertain positive_control}"
export FOLDSEEK_EXTRA_ARGS="${FOLDSEEK_EXTRA_ARGS:--a --alignment-type 1 -e 10 --max-seqs 2000 --threads ${SLURM_CPUS_PER_TASK:-16}}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "QUERY_SETS=${QUERY_SETS} OUT_DIR=${OUT_DIR}"

./structure_prediction/run_foldseek_search.sh

echo "Job finished: $(date)"
