#!/bin/bash
# One-time `mmseqs createdb --createdb-mode 2` on GOPC, building the
# GPU-searchable target database that cluster/run_gopc_search.sbatch's
# --gpu search needs. Separate from cluster/run_gopc_build_db.sbatch (the
# plain/CPU build) and its own separate output path (target_db_gpu/) --
# the two formats are incompatible, so this project keeps them as
# distinct scripts/paths rather than one script with a mode switch, to
# avoid ever being unsure which format a given database is in.
#
# Runs on gpu-a100 (mirroring ../fair_ocean_agent/cluster/run_extraction.sbatch's
# own GPU job), even though this build step itself doesn't strictly need
# CUDA to WRITE the GPU-compatible format (confirmed live: it completes
# fine with no GPU present) -- we have no way to confirm the resulting
# database is actually fully correct for a later --gpu search without a
# real GPU to test against, and getting that wrong after a ~184GB rebuild
# would be far more costly than running this one step on a GPU node we
# didn't strictly need. When in doubt, match the resource type the
# consuming step (the actual --gpu search) needs.
#
# --mem/--cpus-per-task below: `sinfo -p gpu-a100 -o "%P %m %c %l %G"`
# reports each gpu-a100 node as ~1031735MB (~1007GiB) RAM, 128 CPUs, and
# 8 A100 GPUs (gpu:a100:8) -- a big node shared across up to 8 concurrent
# single-GPU jobs. Since we only request --gres=gpu:1, --mem/--cpus-per-task
# below are set to a considerate ~1/8 share (1031735/8 ~ 129000MB,
# 128/8 = 16 CPUs) rather than grabbing disproportionately more of a
# shared node than our 1-GPU request implies -- there's headroom to raise
# --mem further (up to the full ~1TB) if this specific job turns out to
# need it and the cluster isn't busy, but 1/8 is the polite default.
# Confirmed live (see cluster/run_gopc_build_db.sbatch's own header):
# mmseqs REFUSES --no-shuffle together with --createdb-mode 2 and
# force-re-enables shuffle, so this needs MORE memory than the plain
# build's 64G, not the same -- 120G below is a considered guess within
# the fair-share ceiling above, not a validated number.
#
# Don't forget --account on submission if your cluster requires one, e.g.:
#   sbatch --account=191001-364393 cluster/run_gopc_build_db_gpu.sbatch
#
# --time=48:00:00: confirmed live that GPU jobs on this account are capped
# at 48h regardless of gpu-a100's own nominal "infinite" TIMELIMIT (a
# QOS-level restriction tied to the account/allocation, not the raw
# partition config sinfo reports).
#SBATCH --job-name=pha_gopc_build_db_gpu
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gopc_build_db_gpu_%j.out
#SBATCH --error=logs/gopc_build_db_gpu_%j.err

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

# Uses the plain CPU mmseqs build on purpose (this step doesn't need the
# GPU-enabled binary -- --createdb-mode 2 is a storage-format choice, and
# only the actual search step needs mmseqs-linux-gpu).
if ! command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs not found on PATH -- run ./cluster/install_mmseqs2.sh first." >&2
  exit 2
fi
echo "mmseqs: $(mmseqs version)"
echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi || echo "nvidia-smi not available/failed -- is a GPU actually attached to this job? (--gres=gpu:1 above)" >&2

GOPC_FAA="${GOPC_FAA:-PHA_bioprospecting/databases/GOPC/raw/GOPC.geneset.pep.fa.gz}"
if [ ! -f "${GOPC_FAA}" ]; then
  echo "${GOPC_FAA} not found -- download it first:" >&2
  echo "  sbatch --account=191001-364393 --export=ALL,DATASET=gopc,TARGET=geneset cluster/run_download_databases.sbatch" >&2
  exit 2
fi
echo "GOPC_FAA=${GOPC_FAA}"
df -h PHA_bioprospecting/gopc_search 2>/dev/null || df -h .

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/gopc_search/target_db_gpu/gopc_db}"
echo "TARGET_DB=${TARGET_DB}"
echo "NOTE: mmseqs forces --shuffle 1 whenever --createdb-mode 2 is set (confirmed live -- it overrides"
echo "--no-shuffle with a warning: 'Shuffle database cannot be turned off for --createdb-mode 2')."

pha-reference gopc-build-db "${GOPC_FAA}" --target-db "${TARGET_DB}" --threads "${SLURM_CPUS_PER_TASK:-8}" --gpu-compatible

echo "Job finished: $(date)"
echo "GPU-compatible target database ready: ${TARGET_DB}"
echo "Next: cluster/run_gopc_search.sbatch (uses this path by default)"
