#!/bin/bash
# Run TemStaPro on all phaC target sequences as a SLURM array.
#
# First, on service/login node:
#   bash cluster/temstapro/install_temstapro_service.sh /path/to/TemStaPro
#   python3 figures/scripts/export_phac_temstapro_inputs.py
#
# Then submit after checking the chunk count:
#   ls ../fair_ocean_agent/temstapro/chunks/*.faa | wc -l
#   sbatch --array=0-128 --export=ALL,TEMSTAPRO_DIR=/path/to/TemStaPro,TEMSTAPRO_ENV=temstapro_env cluster/temstapro/run_phac_temstapro_array.sbatch
#
# GPU is strongly recommended: TemStaPro README reports ~10 min / 1000 avg
# 1137-aa proteins on an RTX 2080 Ti vs ~10 h on laptop CPU.
#SBATCH --job-name=phac_temstapro
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/phac_temstapro_%A_%a.out
#SBATCH --error=logs/phac_temstapro_%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Run with sbatch from PHA_Ocean_Atlas repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs

: "${TEMSTAPRO_DIR:?Set TEMSTAPRO_DIR to the TemStaPro checkout path.}"
TEMSTAPRO_ENV="${TEMSTAPRO_ENV:-temstapro_env}"

set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${TEMSTAPRO_ENV}"
set -u

BASE="../fair_ocean_agent/temstapro"
CHUNK="${BASE}/chunks/phaC_temstapro_chunk_${SLURM_ARRAY_TASK_ID}.faa"
if [ ! -f "${CHUNK}" ]; then
  CHUNK="$(printf "%s/chunks/phaC_temstapro_chunk_%03d.faa" "${BASE}" "${SLURM_ARRAY_TASK_ID}")"
fi
if [ ! -f "${CHUNK}" ]; then
  echo "Missing chunk for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

OUTDIR="${BASE}/temstapro_outputs"
EMBDIR="${BASE}/embeddings/chunk_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${OUTDIR}" "${EMBDIR}"

OUT="${OUTDIR}/phaC_temstapro_chunk_${SLURM_ARRAY_TASK_ID}.tsv"
echo "started $(date)"
echo "host $(hostname)"
echo "chunk ${CHUNK}"
echo "output ${OUT}"

"${TEMSTAPRO_DIR}/temstapro" \
  --input-fasta "${CHUNK}" \
  --embeddings-dir "${EMBDIR}" \
  --PT-directory "${TEMSTAPRO_DIR}/ProtTrans" \
  --temstapro-directory "${TEMSTAPRO_DIR}" \
  --more-thresholds \
  --mean-output "${OUT}"

echo "finished $(date)"
