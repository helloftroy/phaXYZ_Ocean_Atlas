#!/bin/bash
# Adds real numeric sample depth (meters) to an already-built
# <family>_unique_targets_with_metadata.tsv, by fetching each row's
# sample's underlying NCBI BioSample record -- OMDB's own API
# (omdb-enrich-metadata) has NO depth field at all, confirmed by reading
# their column definitions; depth only exists (for the samples that
# report it) in NCBI's own BioSample metadata. See
# phaatlas/pipeline/ncbi_depth.py's module docstring for the full story,
# including why coverage is inherently partial (not every sample has a
# real NCBI accession; not every accession that does reports depth).
# Needs internet, no GPU -- runs on `service`.
#
# Prerequisite: omdb-enrich-metadata must already have been run for this
# family (produces the *_unique_targets_with_metadata.tsv this reads):
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_enrich_metadata.sbatch
#
# Then:
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_omdb_enrich_depth.sbatch
#
# Follow with the depth-heatmap command directly (fast, local, no
# scheduling needed):
#   pha-reference depth-heatmap --family phaC
#SBATCH --job-name=pha_omdb_enrich_depth
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_enrich_depth_%j.out
#SBATCH --error=logs/omdb_enrich_depth_%j.err

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

: "${FAMILY:?Set FAMILY to a single family_id or comma-separated list.}"
echo "FAMILY=${FAMILY}"

pha-reference omdb-enrich-depth --family "${FAMILY}"

echo "Job finished: $(date)"
echo "Results: PHA_bioprospecting/omdb_search/results/<family>_unique_targets_with_metadata_depth.tsv"
