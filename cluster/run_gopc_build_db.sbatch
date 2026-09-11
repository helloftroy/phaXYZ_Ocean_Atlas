#!/bin/bash
# One-time `mmseqs createdb` on GOPC, building the searchable target
# database that cluster/run_gopc_search.sbatch queries against. Run this
# AFTER downloading GOPC (cluster/run_download_databases.sbatch,
# DATASET=gopc TARGET=geneset) and only once -- the resulting target_db/
# is reused by every family's search afterward.
#
# Heavy I/O over the full ~184GB input, not internet-bound (mmseqs reads
# the local .gz directly, no decompression step needed first) -- give it
# real disk headroom: confirm free scratch space with `df -h
# PHA_bioprospecting/gopc_search` before submitting (the target database
# itself is typically comparable in size to the input, sometimes larger
# once indexed).
#SBATCH --job-name=pha_gopc_build_db
#SBATCH --partition=service
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gopc_build_db_%j.out
#SBATCH --error=logs/gopc_build_db_%j.err

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

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs not found on PATH -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi
echo "mmseqs: $(mmseqs version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"

GOPC_FAA="${GOPC_FAA:-PHA_bioprospecting/databases/GOPC/raw/GOPC.geneset.pep.fa.gz}"
if [ ! -f "${GOPC_FAA}" ]; then
  echo "${GOPC_FAA} not found -- download it first:" >&2
  echo "  sbatch --export=ALL,DATASET=gopc,TARGET=geneset cluster/run_download_databases.sbatch" >&2
  exit 2
fi
echo "GOPC_FAA=${GOPC_FAA}"
df -h PHA_bioprospecting/gopc_search 2>/dev/null || df -h .

pha-reference gopc-build-db "${GOPC_FAA}" --threads "${SLURM_CPUS_PER_TASK:-8}"

echo "Job finished: $(date)"
echo "Target database ready: PHA_bioprospecting/gopc_search/target_db/gopc_db"
echo "Next: cluster/run_gopc_search.sbatch"
