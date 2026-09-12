#!/bin/bash
# One-time `mmseqs createdb` on GOPC, building the PLAIN (CPU-searchable)
# target database. For a GPU-searchable database, use
# cluster/run_gopc_build_db_gpu.sbatch instead -- the two formats are
# incompatible (one built with --gpu-compatible, this one without) and
# this project keeps them as separate scripts/paths on purpose so it's
# never ambiguous which one a given `gopc-search` run needs.
#
# Run this AFTER downloading GOPC (cluster/run_download_databases.sbatch,
# DATASET=gopc TARGET=geneset) and only once -- the resulting target_db/
# is reused by every family's CPU search afterward.
#
# Heavy I/O over the full ~184GB input, not internet-bound (mmseqs reads
# the local .gz directly, no decompression step needed first) -- give it
# real disk headroom: confirm free scratch space with `df -h
# PHA_bioprospecting/gopc_search` before submitting (the target database
# itself is typically comparable in size to the input, sometimes larger
# once indexed).
#
# Don't forget --account on submission if your cluster requires one for
# this partition/QOS (confirmed live: omitting it can silently route the
# job through a different default QOS with its own, possibly tighter,
# memory enforcement -- a job that OOM-kills without --account may run
# fine at the SAME --mem once it's submitted correctly), e.g.:
#   sbatch --account=191001-364393 cluster/run_gopc_build_db.sbatch
#
# SHUFFLE defaults to 0 (`pha-reference gopc-build-db --no-shuffle`, the
# phaatlas default) -- mmseqs createdb's own --shuffle default (on)
# reorders the whole input up front, which costs extra memory at GOPC's
# real scale (hundreds of millions of sequences) for no benefit to this
# one-time build step; --no-shuffle avoids that cost. Re-enable with
# SHUFFLE=1 if a later search's target-split load-balancing seems to
# need it and you have memory to spare.
#
# This job already succeeds at --mem=64G below (with --account and
# --no-shuffle) -- for reference, `sinfo -p service -o "%P %m %c %l"`
# reports this partition's real ceiling as ~128356MB (~125GiB) per node,
# so there's room to raise --mem here too if a future run (e.g. with
# SHUFFLE=1) needs it.
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
  echo "  sbatch --account=191001-364393 --export=ALL,DATASET=gopc,TARGET=geneset cluster/run_download_databases.sbatch" >&2
  exit 2
fi
echo "GOPC_FAA=${GOPC_FAA}"
df -h PHA_bioprospecting/gopc_search 2>/dev/null || df -h .

SHUFFLE_FLAG="--no-shuffle"
if [ "${SHUFFLE:-0}" = "1" ]; then
  SHUFFLE_FLAG="--shuffle"
  echo "SHUFFLE=1 -- re-enabling mmseqs createdb's own shuffle (better target-split load-balancing" \
       "later, at some extra memory cost) -- raise --mem above if this runs low on memory."
fi

TARGET_DB="${TARGET_DB:-PHA_bioprospecting/gopc_search/target_db/gopc_db}"
echo "TARGET_DB=${TARGET_DB}"

pha-reference gopc-build-db "${GOPC_FAA}" --target-db "${TARGET_DB}" --threads "${SLURM_CPUS_PER_TASK:-8}" "${SHUFFLE_FLAG}"

echo "Job finished: $(date)"
echo "Target database ready: ${TARGET_DB}"
echo "Next (CPU search): cluster/run_gopc_search.sbatch --split-memory-limit fallback,"
echo "or build a GPU-compatible database with cluster/run_gopc_build_db_gpu.sbatch instead."
