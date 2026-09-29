#!/bin/bash
# Direct KOfam marker scan for the 5 CBB-positive Thioglobus/SUP05 proteomes.
#
# This bypasses KOfamScan's Ruby/GNU-parallel wrapper and runs hmmsearch
# directly against selected KOfam profiles, using ko_list thresholds.
#
# Submit from PHA_Ocean_Atlas root:
#   sbatch --account=191001-364393 \
#     --export=ALL,KOFAM_ENV=/scratch/morrill/users/hmp278/conda_envs/faire-agent,KOFAM_DB=/scratch/morrill/users/hmp278/db/kofam \
#     cluster/run_thioglobus5_kofam_marker_direct.sbatch
#SBATCH --job-name=thioglobus5_kofam_direct
#SBATCH --partition=service
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/thioglobus5_kofam_direct_%j.out
#SBATCH --error=logs/thioglobus5_kofam_direct_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Run with sbatch from PHA_Ocean_Atlas repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

if [ -n "${KOFAM_ENV:-}" ]; then
  if [ ! -d "${KOFAM_ENV}/bin" ]; then
    echo "KOFAM_ENV was set, but ${KOFAM_ENV}/bin does not exist." >&2
    exit 2
  fi
  export PATH="${KOFAM_ENV}/bin:${PATH}"
fi

: "${KOFAM_DB:?Set KOFAM_DB to the directory containing ko_list and profiles/.}"
if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "hmmsearch not found. Set KOFAM_ENV to an env containing hmmer." >&2
  exit 2
fi

echo "started $(date)"
echo "host $(hostname)"
echo "hmmsearch=$(command -v hmmsearch)"
echo "KOFAM_DB=${KOFAM_DB}"

python3 figures/scripts/run_thioglobus5_kofam_marker_direct.py \
  --kofam-db "${KOFAM_DB}" \
  --hmmsearch-bin "$(command -v hmmsearch)"

echo "finished $(date)"
echo "Pull back:"
echo "  figures/thioglobus5_kofam_marker_summary.tsv"
echo "  figures/thioglobus5_kofam_marker_hits.tsv"
echo "  thioglobus5_reannotation/kofam_marker_direct/"
