#!/bin/bash
# Joins a family's OMDB search results (<family>_unique_targets.tsv) to
# real genome/sample metadata -- GTDB taxonomy per genome, lat/lon and
# ecosystem context per sample -- entirely via OMDB's own data, in two
# steps: (1) stream OMDBv2.0_AA_G_NR100.cluster.tsv.gz once to resolve
# each target_id (an NR100 cluster representative) to its member
# genome(s), (2) batch-query omdb.microbiomics.io's own public JSON API
# (the same one their web genome browser uses -- confirmed live, no auth
# needed) for taxonomy and location. See
# phaatlas/pipeline/omdb_metadata.py's module docstring for the full
# story on both sources, including why OMDB has no raw numeric depth
# field. Needs internet -- runs on `service`, no GPU/mmseqs involved at
# all.
#
# One-time prerequisite: download the cluster file (~4.2GB, not needed
# for the search stage itself, only this one):
#   sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100-clusters \
#     cluster/run_download_databases.sbatch
#
# Start with a SMALL family first, same philosophy as every other stage
# in this pipeline -- phaQ (~5K unique_targets.tsv rows) is the smallest;
# phaB/phaA (500K+ rows) will need many more genome/sample API batches
# and correspondingly longer:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaQ cluster/run_omdb_enrich_metadata.sbatch
# Then a bigger one once you know the timing/coverage looks right:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_enrich_metadata.sbatch
#
# --time=48:00:00: same generous default used for every script in this
# project now, not a measured number for this specific job -- streaming
# the real 249.5M-row cluster file end to end has not been timed yet
# (only tested live against a small real slice of it), but is expected to
# take well under 48h regardless of family size. If a run is killed by
# the time limit, it simply restarts the whole streaming pass on
# resubmit -- there is no partial resume here (unlike the search stage),
# since a single pass over the cluster file is expected to be fast
# relative to search in the first place.
#
# --mem: the per-target_id result dict is small (a handful of strings per
# match), so memory scales with how many target_ids matched, not with the
# cluster file size itself -- 16G is comfortable headroom even for phaB's
# 552K rows.
#SBATCH --job-name=pha_omdb_enrich_metadata
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_enrich_metadata_%j.out
#SBATCH --error=logs/omdb_enrich_metadata_%j.err

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

echo "Job started: $(date)"
echo "Host: $(hostname)"

: "${FAMILY:?Set FAMILY to a single family_id or comma-separated list. Not all, deliberately -- see this script header.}"
CLUSTER_TSV="${CLUSTER_TSV:-PHA_bioprospecting/databases/OMDBv2/OMDBv2.0_AA_G_NR100.cluster.tsv.gz}"
MAX_GENOMES_PER_TARGET="${MAX_GENOMES_PER_TARGET:-5}"
SLEEP_SECONDS="${SLEEP_SECONDS:-0.3}"
echo "FAMILY=${FAMILY} CLUSTER_TSV=${CLUSTER_TSV} MAX_GENOMES_PER_TARGET=${MAX_GENOMES_PER_TARGET}"

if [ ! -f "${CLUSTER_TSV}" ]; then
  echo "${CLUSTER_TSV} not found -- download it first:" >&2
  echo "  sbatch --account=191001-364393 --export=ALL,DATASET=omdb,TARGET=nr100-clusters cluster/run_download_databases.sbatch" >&2
  exit 2
fi

pha-reference omdb-enrich-metadata \
  --family "${FAMILY}" \
  --cluster-tsv "${CLUSTER_TSV}" \
  --max-genomes-per-target "${MAX_GENOMES_PER_TARGET}" \
  --sleep-seconds "${SLEEP_SECONDS}"

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/omdb_search/results/<family>_unique_targets_with_metadata.tsv"
