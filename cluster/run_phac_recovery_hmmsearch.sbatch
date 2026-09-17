#!/bin/bash
# Step 3 of 3. Searches phac_recovery/neighborhoods.faa (step 2's output)
# with both phaC HMM profiles at a relaxed --domE 10 threshold (see
# phac_recovery/hmm/ -- Pfam PF07167 N-terminal domain, and phaC_custom.hmm
# built from this repo's own 822-sequence curated phaC diversity). CPU
# only, no internet needed.
#
# One-time prerequisite: hmmer isn't in the pha-reference conda env by
# default -- install it once before the first run of this script:
#   conda install -n pha-reference -c bioconda -c conda-forge hmmer
#
#   sbatch --account=191001-364393 cluster/run_phac_recovery_hmmsearch.sbatch
# Tighten/loosen the threshold if 10 is producing too many/few candidates:
#   sbatch --account=191001-364393 --export=ALL,DOM_E=1 cluster/run_phac_recovery_hmmsearch.sbatch
#SBATCH --job-name=phac_recovery_hmmsearch
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/phac_recovery_hmmsearch_%j.out
#SBATCH --error=logs/phac_recovery_hmmsearch_%j.err

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

if ! command -v hmmsearch >/dev/null 2>&1; then
  echo "hmmsearch not on PATH -- run: conda install -n pha-reference -c bioconda -c conda-forge hmmer" >&2
  exit 2
fi

DOM_E="${DOM_E:-10}"
NEIGHBORHOODS="phac_recovery/neighborhoods.faa"
OUT_HITS="phac_recovery/candidate_phaC_hits.tsv"

if [ ! -f "${NEIGHBORHOODS}" ]; then
  echo "${NEIGHBORHOODS} not found -- run cluster/run_phac_recovery_extract_neighborhoods.sbatch (step 2) first." >&2
  exit 2
fi

pha-reference phac-recovery-hmmsearch \
  --neighborhoods-fasta "${NEIGHBORHOODS}" \
  --out-hits "${OUT_HITS}" \
  --dom-e "${DOM_E}"

echo "Job finished: $(date)"
echo "Candidate phaC hits: ${OUT_HITS}"
echo "scp this back: scp <cluster>:<repo>/${OUT_HITS} fair_ocean_agent/phac_recovery/"
