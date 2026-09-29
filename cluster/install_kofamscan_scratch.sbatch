#!/bin/bash
# Run the scratch-safe KOfamScan/KOfam database installer on a service node.
#
# Submit from the PHA_Ocean_Atlas repo root:
#   sbatch --account=191001-364393 cluster/install_kofamscan_scratch.sbatch
#
# Optional overrides:
#   sbatch --account=191001-364393 \
#     --export=ALL,CONDA_ENV_PREFIX=/scratch/morrill/users/hmp278/conda_envs/faire-agent,KOFAM_DB=/scratch/morrill/users/hmp278/db/kofam \
#     cluster/install_kofamscan_scratch.sbatch
#SBATCH --job-name=kofam_install
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/kofam_install_%j.out
#SBATCH --error=logs/kofam_install_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Run with sbatch from PHA_Ocean_Atlas repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

echo "started $(date)"
echo "host $(hostname)"
echo "submit dir ${SLURM_SUBMIT_DIR}"

bash cluster/install_kofamscan_scratch.sh

echo "finished $(date)"
echo "Next:"
echo "  sbatch --account=191001-364393 --export=ALL,KOFAM_ENV=/scratch/morrill/users/hmp278/conda_envs/faire-agent,KOFAM_DB=/scratch/morrill/users/hmp278/db/kofam cluster/run_thioglobus5_kofamscan.sbatch"
