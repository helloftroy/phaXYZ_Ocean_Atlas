#!/bin/bash
# One-time `mmseqs createdb --createdb-mode 2` on OMDBv2.0_AA_G_NR100,
# building the GPU-searchable target database that
# cluster/run_omdb_search.sbatch's --gpu search needs. Mirrors
# cluster/run_gopc_build_db_gpu.sbatch exactly, just pointed at OMDB.
#
# OMDBv2.0_AA_G_NR100 is ~10x fewer sequences than GOPC (249.5M vs.
# ~2.46B) -- this build and the search itself should need noticeably less
# time/memory than GOPC's did, but that's an expectation, not something
# confirmed live yet. --mem/--cpus-per-task below are left at the SAME
# proven-working values as GOPC's build (see that script's own header for
# the sinfo-derived reasoning) rather than guessed downward -- safer to
# over-provision a shared resource once more than to repeat the OOM/
# scheduling churn already hit twice on GOPC. Tune down once real OMDB
# numbers are in hand.
#
# Don't forget --account on submission if your cluster requires one, e.g.:
#   sbatch --account=191001-364393 cluster/run_omdb_build_db_gpu.sbatch
#
# --time=48:00:00: same account-level GPU cap confirmed live for the GOPC
# jobs applies here too (a QOS restriction, not gpu-a100's own nominal
# "infinite" TIMELIMIT).
#SBATCH --job-name=pha_omdb_build_db_gpu
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_build_db_gpu_%j.out
#SBATCH --error=logs/omdb_build_db_gpu_%j.err

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

# Uses the plain CPU mmseqs build on purpose -- --createdb-mode 2 is a
# storage-format choice, not a CUDA requirement to write; only the actual
# search step needs mmseqs-linux-gpu.
if ! command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs not found on PATH -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi
echo "mmseqs: $(mmseqs version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi || echo "nvidia-smi not available/failed -- is a GPU actually attached to this job? (--gres=gpu:1 above)" >&2

OMDB_FAA="${OMDB_FAA:-PHA_bioprospecting/databases/OMDBv2/OMDBv2.0_AA_G_NR100.faa.gz}"
if [ ! -f "${OMDB_FAA}" ]; then
  echo "${OMDB_FAA} not found -- download it first:" >&2
  echo "  sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100 cluster/run_download_databases.sbatch" >&2
  exit 2
fi
echo "OMDB_FAA=${OMDB_FAA}"
df -h PHA_bioprospecting/omdb_search 2>/dev/null || df -h .

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/omdb_search/target_db_gpu/omdb_db}"
echo "TARGET_DB=${TARGET_DB}"
echo "NOTE: mmseqs forces --shuffle 1 whenever --createdb-mode 2 is set (confirmed live on GOPC -- it overrides"
echo "--no-shuffle with a warning: 'Shuffle database cannot be turned off for --createdb-mode 2')."

pha-reference omdb-build-db "${OMDB_FAA}" --target-db "${TARGET_DB}" --threads "${SLURM_CPUS_PER_TASK:-8}" --gpu-compatible

echo "Job finished: $(date)"
echo "GPU-compatible target database ready: ${TARGET_DB}"
echo "Next: cluster/run_omdb_search.sbatch (uses this path by default)"
