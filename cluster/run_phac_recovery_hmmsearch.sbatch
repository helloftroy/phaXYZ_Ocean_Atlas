#!/bin/bash
# Step 3 of 3. Searches phac_recovery/neighborhoods.faa (step 2's output)
# with both phaC HMM profiles at a relaxed --domE 10 threshold (see
# phac_recovery/hmm/ -- Pfam PF07167 N-terminal domain, and phaC_custom.hmm
# built from this repo's own 822-sequence curated phaC diversity). CPU
# only, no internet needed.
#
# One-time prerequisite: hmmer needs to be installed somewhere. Installing
# it INTO an existing env (e.g. faire-agent) can fail with "Disk quota
# exceeded" writing a .condatmp dir even after redirecting CONDA_PKGS_DIRS
# to scratch -- confirmed live -- because that existing env's own folder
# still lives under $HOME/miniforge3/envs/, and conda needs to write inside
# the env's own prefix during install regardless of where the package cache
# points. Simplest fix: don't touch that env at all -- create a small,
# standalone env for just hmmer, entirely on scratch:
#   export CONDA_REGISTER_ENVS=false   # skip even the tiny ~/.conda/environments.txt write
#   conda create -p /scratch/$USER/conda_envs/hmmer_env -c bioconda -c conda-forge hmmer -y
# Then point this script at it (it's not on PATH, so pass it explicitly):
#   sbatch --account=191001-364393 --export=ALL,HMMER_BIN_DIR=/scratch/$USER/conda_envs/hmmer_env/bin cluster/run_phac_recovery_hmmsearch.sbatch
#
# Tighten/loosen the threshold if 10 is producing too many/few candidates:
#   sbatch --account=191001-364393 --export=ALL,HMMER_BIN_DIR=...,DOM_E=1 cluster/run_phac_recovery_hmmsearch.sbatch
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

HMMSEARCH_BIN="hmmsearch"
if [ -n "${HMMER_BIN_DIR:-}" ]; then
  HMMSEARCH_BIN="${HMMER_BIN_DIR}/hmmsearch"
fi
if ! command -v "${HMMSEARCH_BIN}" >/dev/null 2>&1; then
  echo "${HMMSEARCH_BIN} not found/executable -- install hmmer first (see this script's header comment for the scratch-env fix) and pass its bin dir via HMMER_BIN_DIR." >&2
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
  --dom-e "${DOM_E}" \
  --hmmsearch-bin "${HMMSEARCH_BIN}"

echo "Job finished: $(date)"
echo "Candidate phaC hits: ${OUT_HITS}"
echo "scp this back: scp <cluster>:<repo>/${OUT_HITS} fair_ocean_agent/phac_recovery/"
