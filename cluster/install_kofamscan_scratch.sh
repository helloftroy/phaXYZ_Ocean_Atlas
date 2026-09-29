#!/bin/bash
# Install KOfamScan and download the KOfam database without writing to $HOME.
#
# Intended cluster use:
#   bash cluster/install_kofamscan_scratch.sh
#
# Defaults are set for hmp278's Morrill scratch layout, but can be overridden:
#   CONDA_ENV_PREFIX=/scratch/morrill/users/hmp278/conda_envs/faire-agent \
#   KOFAM_DB=/scratch/morrill/users/hmp278/db/kofam \
#   bash cluster/install_kofamscan_scratch.sh

set -euo pipefail

CONDA_ENV_PREFIX="${CONDA_ENV_PREFIX:-/scratch/morrill/users/hmp278/conda_envs/faire-agent}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/scratch/morrill/users/hmp278}"
CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${SCRATCH_ROOT}/conda_pkgs}"
CONDA_ENVS_PATH="${CONDA_ENVS_PATH:-${SCRATCH_ROOT}/conda_envs}"
KOFAM_DB="${KOFAM_DB:-${SCRATCH_ROOT}/db/kofam}"

export CONDA_PKGS_DIRS
export CONDA_ENVS_PATH
export CONDA_REGISTER_ENVS=false

mkdir -p "${CONDA_PKGS_DIRS}" "${CONDA_ENVS_PATH}" "${KOFAM_DB}"

echo "CONDA_ENV_PREFIX=${CONDA_ENV_PREFIX}"
echo "CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS}"
echo "KOFAM_DB=${KOFAM_DB}"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found on PATH. Load/initialize conda first, then rerun." >&2
  exit 2
fi

if [ ! -d "${CONDA_ENV_PREFIX}" ]; then
  echo "Conda env prefix does not exist yet; creating it in scratch."
  conda create -p "${CONDA_ENV_PREFIX}" -c bioconda -c conda-forge kofamscan -y
else
  echo "Installing/updating KOfamScan in existing scratch env."
  conda install -p "${CONDA_ENV_PREFIX}" -c bioconda -c conda-forge kofamscan -y
fi

download() {
  local url="$1"
  local out="$2"
  if [ -s "${out}" ]; then
    echo "Already present: ${out}"
    return
  fi
  if command -v wget >/dev/null 2>&1; then
    wget -O "${out}" "${url}"
  else
    curl -L --fail --show-error -o "${out}" "${url}"
  fi
}

cd "${KOFAM_DB}"
download "https://www.genome.jp/ftp/db/kofam/ko_list.gz" "ko_list.gz"
download "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz" "profiles.tar.gz"

if [ ! -f ko_list ]; then
  gunzip -kf ko_list.gz
fi
if [ ! -d profiles ]; then
  tar -xzf profiles.tar.gz
fi

echo "Done."
echo "Use with:"
echo "  sbatch --account=191001-364393 --export=ALL,KOFAM_ENV=${CONDA_ENV_PREFIX},KOFAM_DB=${KOFAM_DB} cluster/run_thioglobus5_kofamscan.sbatch"
