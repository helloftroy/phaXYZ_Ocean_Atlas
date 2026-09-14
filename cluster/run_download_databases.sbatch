#!/bin/bash
# Downloads GOPC (~184GB) and/or OMDBv2 (~49GB + ~4.2GB + ~7.5MB) via
# PHA_bioprospecting/scripts/download_{gopc,omdb}.sh -- resumable (curl
# -C -), checksum-verified where a checksum is published, manifest written
# either way. Needs internet (service/HTTPS-capable partition, no GPU).
#
# Confirm free scratch space before submitting -- `df -h` on wherever
# PHA_bioprospecting/databases/ actually lives (check with `df -h
# PHA_bioprospecting/databases` from the repo root) needs on the order of
# 240GB free for everything, or just the relevant subset for TARGET below.
#
# TARGET selects what to fetch (see the two scripts' own --help-equivalent
# usage comments for the full option list). Add --account=<your account>
# if your cluster requires one for this partition/QOS (confirmed live:
# omitting it can silently route the job through a different default QOS
# with its own, possibly tighter, resource enforcement):
#   sbatch --account=191001-364393 --export=ALL,DATASET=gopc,TARGET=md5              cluster/run_download_databases.sbatch
#   sbatch --account=191001-364393 --export=ALL,DATASET=gopc,TARGET=geneset           cluster/run_download_databases.sbatch   # ~184GB
#   sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=data               cluster/run_download_databases.sbatch
#   sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100               cluster/run_download_databases.sbatch  # ~49GB
#   sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100-clusters      cluster/run_download_databases.sbatch  # ~4.2GB
#   sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=all                 cluster/run_download_databases.sbatch
#
# --time=0 (sbatch's own syntax for "no limit") since these are genuinely
# large transfers over a real network -- `sinfo -p service -o "%P %m %c
# %l"` confirmed this partition's own TIMELIMIT is "infinite", so there's
# no cluster-side reason to cap this ourselves. If it's still interrupted
# somehow (node reboot, etc.), just resubmit the exact same command: curl
# -C - resumes from wherever the partial file stopped rather than restarting.
#SBATCH --job-name=pha_download_databases
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0
#SBATCH --output=logs/download_databases_%j.out
#SBATCH --error=logs/download_databases_%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/)
mkdir -p logs

: "${DATASET:?Set DATASET=gopc or DATASET=omdb (see this script's header comment for examples).}"
: "${TARGET:?Set TARGET (see this script's header comment for the valid values per DATASET).}"

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "DATASET=${DATASET} TARGET=${TARGET}"
df -h PHA_bioprospecting/databases 2>/dev/null || df -h .

case "${DATASET}" in
  gopc) ./PHA_bioprospecting/scripts/download_gopc.sh "${TARGET}" ;;
  omdb) ./PHA_bioprospecting/scripts/download_omdb.sh "${TARGET}" ;;
  *) echo "Unknown DATASET='${DATASET}'. Use gopc or omdb." >&2; exit 2 ;;
esac

echo "Job finished: $(date)"
