#!/bin/bash
# Builds a 2D sequence-space embedding of a family's extracted+clustered
# proteins (PCA always, UMAP if installed) -- "each point is a protein",
# colorable by sequence cluster, depth zone, taxon, or pathway
# architecture. CPU-only, no GPU needed -- runs on `service`.
#
# Exists specifically to run the FULL 3-mer (8,000-dim) feature matrix at
# k=3 (pipeline/sequence_embedding.py's own documented default) plus a
# real UMAP embedding -- both were skipped in a local, 16GB-RAM
# exploratory pass (used k=2/400-dim there, and umap-learn was not
# installed) since the raw feature matrix alone is ~8GB at k=3 for
# ~128K sequences, before PCA/UMAP's own working-memory overhead. The
# cluster has far more headroom for both.
#
# One-time prerequisite: install the optional umap-learn extra (heavier,
# pulls in numba -- not part of the base install):
#   conda activate <your env>   # or: source cluster/env_activate.sh
#   pip install -e '.[umap]'
#
# Prerequisites (from cluster/run_omdb_extract_cluster_sequences.sbatch /
# run_omdb_cluster_sequences.sbatch): the extracted FASTA and a cluster
# assignment TSV for this family, at whatever identity threshold you want
# the "cluster" coloring in the output to reflect (e.g. 60%):
#   sbatch --account=191001-364393 --export=ALL,FAMILY=phaC,FASTA=PHA_bioprospecting/omdb_search/results/phaC_cluster_sequences.faa,CLUSTER_TSV=PHA_bioprospecting/omdb_search/results/phaC_cluster0.6_cluster.tsv \
#     cluster/run_omdb_sequence_embedding.sbatch
#
# Output is named by k (K defaults to 3 here) so it never collides with
# a different-k run's output -- e.g. a prior local k=2 exploratory pass.
#SBATCH --job-name=pha_omdb_sequence_embedding
#SBATCH --partition=service
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/omdb_sequence_embedding_%j.out
#SBATCH --error=logs/omdb_sequence_embedding_%j.err

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

: "${FAMILY:?Set FAMILY to a single family_id.}"
: "${FASTA:?Set FASTA to the extracted-sequences fasta path.}"
: "${CLUSTER_TSV:?Set CLUSTER_TSV to a cluster assignment tsv path.}"
K="${K:-3}"
RESULTS_DIR="${RESULTS_DIR:-PHA_bioprospecting/omdb_search/results}"
OUT_PATH="${OUT_PATH:-${RESULTS_DIR}/${FAMILY}_sequence_embedding_k${K}.tsv}"
echo "FAMILY=${FAMILY} FASTA=${FASTA} CLUSTER_TSV=${CLUSTER_TSV} K=${K}"

python3 -c "import umap" 2>/dev/null \
  && echo "umap-learn is installed -- will compute a real UMAP embedding alongside PCA." \
  || echo "umap-learn NOT installed -- UMAP columns will be left blank (PCA still computed). See this script's header for the one-time install command."

pha-reference sequence-embedding \
  --family "${FAMILY}" \
  --fasta "${FASTA}" \
  --cluster-tsv "${CLUSTER_TSV}" \
  --results-dir "${RESULTS_DIR}" \
  --k "${K}" \
  --umap \
  --out-path "${OUT_PATH}"

echo "Job finished: $(date)"
echo "Embedding: ${OUT_PATH}"
echo "scp this file back, e.g.:"
echo "  scp <user>@<cluster>:$(pwd)/${OUT_PATH} ."
