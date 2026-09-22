#!/usr/bin/env bash
# Structural comparison step from PHA_CLEAN_RESULTS.md section 6.3: builds
# a Foldseek database from the folded, audited reference structures, then
# searches the uncertain-tier and positive-control structures against it.
# Run AFTER all three cluster/run_esmfold.sbatch array jobs have finished
# (or far enough along to test on) -- needs real .pdb files under
# structure_prediction/esmfold_out/{reference,uncertain,positive_control}/.
#
# Foldseek's structural alignment (not sequence) is the actual "does this
# fold look like a real PhaC catalytic domain" test -- independent of every
# sequence-based check already run (HMM, mmseqs homology, catalytic-triad
# projection, ESMFold's own pLDDT). A query hitting a reference with a high
# TM-score/high LDDT is real structural evidence, regardless of how low its
# sequence identity was.
#
# Usage (after ESMFold has produced at least the reference + one query set):
#   ./structure_prediction/run_foldseek_search.sh
#   # or override any of the defaults below, e.g. to test on just the
#   # positive controls before the full uncertain-tier run finishes:
#   QUERY_SETS=positive_control ./structure_prediction/run_foldseek_search.sh
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root (PHA_Ocean_Atlas/)

FOLDSEEK_BIN="${FOLDSEEK_BIN:-cluster/bin/foldseek}"
ESMFOLD_OUT="${ESMFOLD_OUT:-structure_prediction/esmfold_out}"
OUT_DIR="${OUT_DIR:-structure_prediction/foldseek_out}"
QUERY_SETS="${QUERY_SETS:-uncertain positive_control}"
# TM-align-based alignment (-a) is slower than the default 3Di+AA search but
# gives TM-score directly, the metric this comparison actually needs --
# worth the extra time at this dataset's scale (~7,300 structures total).
FOLDSEEK_EXTRA_ARGS="${FOLDSEEK_EXTRA_ARGS:--a --alignment-type 1}"

if [ ! -x "${FOLDSEEK_BIN}" ]; then
  echo "${FOLDSEEK_BIN} not found/executable -- run ./cluster/install_foldseek.sh first." >&2
  exit 2
fi

REF_DIR="${ESMFOLD_OUT}/reference"
if [ ! -d "${REF_DIR}" ] || [ -z "$(ls -A "${REF_DIR}"/*.pdb 2>/dev/null)" ]; then
  echo "${REF_DIR} has no .pdb files yet -- the reference-set ESMFold array job" >&2
  echo "needs to finish (or at least partially finish) before this can run." >&2
  exit 2
fi

mkdir -p "${OUT_DIR}/tmp"

echo "Building Foldseek reference database from ${REF_DIR} ..."
"${FOLDSEEK_BIN}" createdb "${REF_DIR}" "${OUT_DIR}/reference_db"

for query_set in ${QUERY_SETS}; do
  QUERY_DIR="${ESMFOLD_OUT}/${query_set}"
  if [ ! -d "${QUERY_DIR}" ] || [ -z "$(ls -A "${QUERY_DIR}"/*.pdb 2>/dev/null)" ]; then
    echo "Skipping ${query_set}: ${QUERY_DIR} has no .pdb files yet."
    continue
  fi
  N_QUERIES=$(ls "${QUERY_DIR}"/*.pdb | wc -l)
  echo
  echo "Searching ${query_set} (${N_QUERIES} structures) against the reference database ..."
  "${FOLDSEEK_BIN}" easy-search "${QUERY_DIR}" "${OUT_DIR}/reference_db" \
    "${OUT_DIR}/${query_set}_vs_reference.tsv" "${OUT_DIR}/tmp" \
    ${FOLDSEEK_EXTRA_ARGS} \
    --format-output "query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,prob,qlen,tlen,alnlen"
  echo "wrote ${OUT_DIR}/${query_set}_vs_reference.tsv"
done

echo
echo "Done. Columns in each output TSV: query, target, evalue, bits, alntmscore (TM-score"
echo "normalized by the alignment), qtmscore/ttmscore (normalized by query/target length"
echo "respectively -- ttmscore is usually the more comparable one across queries of very"
echo "different lengths), lddt, prob (Foldseek's own hit-confidence probability), qlen,"
echo "tlen, alnlen. A high ttmscore/lddt against ANY reference is the positive signal to"
echo "look for -- see PHA_CLEAN_RESULTS.md section 6.3 for how to read these against the"
echo "positive-control set before trusting the uncertain-tier results."
