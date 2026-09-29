#!/bin/bash
# Run METABOLIC-G on the 5 CBB-positive Thioglobus/SUP05 genome FASTAs.
#
# This uses genome FASTAs rather than the protein FASTAs because METABOLIC-G's
# standard MAG workflow expects genome bins and handles its own gene calls /
# metabolic module summaries.
#
# Prepare inputs first:
#   python3 figures/scripts/export_thioglobus5_reannotation_inputs.py
#
# Submit from PHA_Ocean_Atlas root:
#   sbatch --account=191001-364393 \
#     --export=ALL,METABOLIC_DIR=/scratch/$USER/METABOLIC \
#     cluster/run_thioglobus5_metabolic.sbatch
#SBATCH --job-name=thioglobus5_metabolic
#SBATCH --partition=service
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/thioglobus5_metabolic_%j.out
#SBATCH --error=logs/thioglobus5_metabolic_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Run with sbatch from PHA_Ocean_Atlas repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs thioglobus5_reannotation/metabolic

: "${METABOLIC_DIR:?Set METABOLIC_DIR to the METABOLIC checkout/install directory.}"
METABOLIC_G="${METABOLIC_G:-${METABOLIC_DIR}/METABOLIC-G.pl}"
if [ ! -f "${METABOLIC_G}" ]; then
  echo "METABOLIC-G.pl not found at ${METABOLIC_G}" >&2
  exit 2
fi

echo "started $(date)"
echo "host $(hostname)"
echo "METABOLIC_G=${METABOLIC_G}"

perl "${METABOLIC_G}" \
  -in-gn thioglobus5_reannotation/genomes \
  -o thioglobus5_reannotation/metabolic \
  -t "${SLURM_CPUS_PER_TASK:-16}"

echo "finished $(date)"
echo "Pull back: thioglobus5_reannotation/metabolic/"
