#!/bin/bash
# Target-side data prep for the mmseqs-based novel-phaC search: downloads
# the FULL gene-call set (not just anchor-relative neighborhoods -- see
# cluster/run_phac_recovery_extract_neighborhoods.sbatch for that older,
# HMM-based recovery attempt) for all 7,826 qualifying_genomes.txt genomes
# (no phaC hit, >5 other PHA-pathway families present), so the mmseqs
# search doesn't rely on the missed synthase sitting near a known anchor
# gene.
#
# Query-side prep (run locally first, no cluster needed):
#   .venv/bin/python phac_recovery/build_novel_phac_query_set.py
#
# NEEDS INTERNET (7,826 individual HTTPS fetches to sunagawalab.ethz.ch) --
# same partition as cluster/run_phac_recovery_extract_neighborhoods.sbatch
# and cluster/run_omdb_enrich_metadata.sbatch, both of which do the same
# kind of per-genome HTTPS fetch. NOTE: cluster/run_foldseek.sbatch found
# that "service" is this codebase's own placeholder name and the real
# CPU-compute partition for this account is `morrill` -- but that finding
# was for a pure-compute, no-internet job. This step's bottleneck is
# internet access, a different axis than CPU vs GPU, and whether `morrill`
# has outbound internet is NOT yet confirmed -- if this job fails to
# resolve sunagawalab.ethz.ch on morrill, fall back to `service` (or
# whatever `sinfo`/a cluster admin confirms is the internet-capable
# partition for this account) rather than assuming morrill works here too.
# Long walltime because of the sheer number of individual requests, not
# because any one is slow -- 47h is the longest this cluster allows (48h
# cap); if it doesn't finish in that window, resubmitting is NOT
# idempotent yet (always starts from genome_urls[0] again, no
# checkpoint/resume).
#
#   sbatch --account=191001-364393 cluster/run_phac_recovery_download_proteomes.sbatch
#SBATCH --job-name=phac_recovery_download_proteomes
#SBATCH --partition=service
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=47:00:00
#SBATCH --output=logs/phac_recovery_download_proteomes_%j.out
#SBATCH --error=logs/phac_recovery_download_proteomes_%j.err

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

MANIFEST="phac_recovery/genome_download_manifest.tsv"
OUT_FASTA="phac_recovery/novel_phac_target_proteomes.faa"
OUT_REPORT="phac_recovery/proteome_download_report.tsv"

if [ ! -f "${MANIFEST}" ]; then
  echo "${MANIFEST} not found -- run phac_recovery/build_local_manifests.py first (it is already committed/done for this repo's current 7,826-genome set)." >&2
  exit 2
fi

pha-reference phac-recovery-download-proteomes \
  --genome-download-manifest "${MANIFEST}" \
  --out-fasta "${OUT_FASTA}" \
  --out-report "${OUT_REPORT}"

echo "Job finished: $(date)"
echo "Next: sbatch cluster/run_novel_phac_mmseqs_search.sbatch"
