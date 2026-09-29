#!/bin/bash
# Run KOfamScan on the 5 CBB-positive Thioglobus/SUP05 protein FASTAs.
#
# First prepare inputs locally or on the cluster repo copy:
#   python3 figures/scripts/export_thioglobus5_reannotation_inputs.py
#
# One-time KOfam setup on cluster scratch, if needed. This avoids $HOME:
#   bash cluster/install_kofamscan_scratch.sh
#
# Submit from PHA_Ocean_Atlas root on a service node:
#   sbatch --account=191001-364393 \
#     --export=ALL,KOFAM_ENV=/scratch/morrill/users/hmp278/conda_envs/faire-agent,KOFAM_DB=/scratch/morrill/users/hmp278/db/kofam \
#     cluster/run_thioglobus5_kofamscan.sbatch
#SBATCH --job-name=thioglobus5_kofam
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/thioglobus5_kofam_%j.out
#SBATCH --error=logs/thioglobus5_kofam_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Run with sbatch from PHA_Ocean_Atlas repo root.}"
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs thioglobus5_reannotation/kofam thioglobus5_reannotation/kofam_tmp

if [ -n "${KOFAM_ENV:-}" ]; then
  if [ ! -d "${KOFAM_ENV}/bin" ]; then
    echo "KOFAM_ENV was set, but ${KOFAM_ENV}/bin does not exist." >&2
    exit 2
  fi
  export PATH="${KOFAM_ENV}/bin:${PATH}"
fi
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
KOFAM_CPUS="${KOFAM_CPUS:-1}"

: "${KOFAM_DB:?Set KOFAM_DB to the directory containing ko_list and profiles/.}"
KOFAM_PROFILES="${KOFAM_PROFILES:-${KOFAM_DB}/profiles}"
KOFAM_KO_LIST="${KOFAM_KO_LIST:-${KOFAM_DB}/ko_list}"

if ! command -v exec_annotation >/dev/null 2>&1; then
  echo "exec_annotation not found. Set KOFAM_ENV or put KOfamScan on PATH." >&2
  exit 2
fi
if [ ! -d "${KOFAM_PROFILES}" ] || [ ! -f "${KOFAM_KO_LIST}" ]; then
  echo "Missing KOfam profiles or ko_list under KOFAM_DB=${KOFAM_DB}" >&2
  exit 2
fi

echo "started $(date)"
echo "host $(hostname)"
echo "KOFAM_PROFILES=${KOFAM_PROFILES}"
echo "KOFAM_KO_LIST=${KOFAM_KO_LIST}"
echo "KOFAM_CPUS=${KOFAM_CPUS}"
echo "exec_annotation=$(command -v exec_annotation)"
echo "hmmsearch=$(command -v hmmsearch || true)"
echo "parallel=$(command -v parallel || true)"
exec_annotation --version || true

for faa in thioglobus5_reannotation/proteins/*.faa; do
  genome="$(basename "${faa}" .faa)"
  out="thioglobus5_reannotation/kofam/${genome}.kofam.detail.tsv"
  tmp="thioglobus5_reannotation/kofam_tmp/${genome}"
  rm -rf "${tmp}"
  mkdir -p "${tmp}"
  echo "KOfamScan ${genome}"
  exec_annotation \
    -p "${KOFAM_PROFILES}" \
    -k "${KOFAM_KO_LIST}" \
    --cpu "${KOFAM_CPUS}" \
    --tmp-dir "${tmp}" \
    -f detail-tsv \
    -o "${out}" \
    "${faa}"
done

python3 figures/scripts/summarize_thioglobus5_kofam.py

echo "finished $(date)"
echo "Pull back:"
echo "  thioglobus5_reannotation/kofam/"
echo "  figures/thioglobus5_kofam_marker_summary.tsv"
echo "  figures/thioglobus5_kofam_marker_hits.tsv"
