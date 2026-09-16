#!/bin/bash
# Builds/refreshes PHA_reference/pha_reference.sqlite for every family in
# config/family_definitions.yaml: BRENDA download + EC-anchored ingest,
# then UniProt Route A/B ingest, then export. No GPU needed -- this is
# pure HTTPS retrieval (BRENDA + UniProt REST), so it belongs on the same
# kind of CPU/HTTPS-capable partition as ../fair_ocean_agent's own
# run_discovery.sbatch, not a GPU one.
#
# Partition name below is copied from that same convention -- confirm the
# real partition name for your account with `sinfo` before submitting.
#
# --time=48:00:00: same generous default used for every script in this
# project now, not a cap this job is expected to need -- `sinfo -p
# service` confirmed this partition's own TIMELIMIT is "infinite" too.
#SBATCH --job-name=pha_reference_build
#SBATCH --partition=service
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pha_reference_build_%j.out
#SBATCH --error=logs/pha_reference_build_%j.err

set -euo pipefail
# SLURM_SUBMIT_DIR, not $(dirname "$0") -- sbatch commonly stages/copies
# the submitted script before running it, so $0 can point at a spool
# directory with no write access (see ../fair_ocean_agent/cluster's own
# sbatch scripts, where this was confirmed live).
: "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset -- this script must be run via sbatch, not directly.}"
cd "${SLURM_SUBMIT_DIR}"   # repo root (PHA_Ocean_Atlas/), i.e. wherever you ran `sbatch cluster/run_pha_reference_build.sbatch` from
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
echo "Working directory: $(pwd)"
python --version

if [ ! -f .env ] && [ -f .env.example ]; then
  echo "No .env found -- BRENDA/UniProt requests will use the placeholder" >&2
  echo "contact email. cp .env.example .env and set PHA_REFERENCE_CONTACT_EMAIL" >&2
  echo "before a large run (see cluster/README.md)." >&2
fi

# FAMILY selects which family_id(s) to run (comma-separated, or 'all' --
# the default). Test with a single family first, e.g. (add
# --account=<your account> if your cluster requires one for this
# partition/QOS):
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC cluster/run_pha_reference_build.sbatch
FAMILY="${FAMILY:-all}"
echo "FAMILY=${FAMILY}"

pha-reference init-db

echo "=== BRENDA: EC-anchored ingest (families with ec_brenda set; download is cached after the first run) ==="
pha-reference ingest-brenda --family "${FAMILY}"

echo "=== UniProt: Route B family-specific search (also runs BRENDA-linked Route A fetches inline) ==="
pha-reference ingest-uniprot --family "${FAMILY}"

echo "=== export ==="
pha-reference export
pha-reference status

echo "Job finished: $(date)"
echo "Results: PHA_reference/pha_reference.sqlite, PHA_reference/exports/"
echo "scp that directory back, e.g.:"
echo "  scp -r <cluster>:$(pwd)/PHA_reference ./"
