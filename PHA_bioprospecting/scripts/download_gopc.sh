#!/usr/bin/env bash
# Downloads GOPC (Global Ocean Protein Catalog), the gene-set FASTA output
# of the Global Ocean Microbiome Catalogue (GOMC) -- Chen, Jia, Sun et al.,
# "Global marine microbial diversity and its potential in bioprospecting",
# Nature 633:371-379 (2024), doi:10.1038/s41586-024-07891-2. Hosted on
# CNGBdb as dataset MDB0000002; the actual files live on CNGB's own public
# FTP-over-HTTPS mirror (no login), confirmed to support HTTP Range
# requests (resumable).
#
# Usage:
#   ./download_gopc.sh md5        # md5.txt only (tiny, safe to run anywhere)
#   ./download_gopc.sh geneset    # GOPC.geneset.pep.fa.gz -- ~184GB. Cluster only.
#   ./download_gopc.sh all        # both
#
# GOPC.geneset.pep.fa.gz is ~184GB (183,959,496,042 bytes at time of
# writing) -- confirm free scratch space before running the `geneset`
# target (df -h on the target filesystem). This script does not decompress
# it; leave it as .gz until a downstream search tool actually needs the
# unpacked form (see PHA_bioprospecting's own README for why).
set -euo pipefail
cd "$(dirname "$0")"
source ./lib_download.sh

BASE_URL="https://ftp.cngb.org/pub/SciRAID/microbiomics/MDB0000002"
RAW_DIR="../databases/GOPC/raw"
MANIFEST="../databases/GOPC/download_manifest.tsv"

TARGET="${1:-}"
if [ -z "${TARGET}" ]; then
  echo "Usage: $0 {md5|geneset|all}" >&2
  exit 2
fi

download_md5() {
  local dest="${RAW_DIR}/md5.txt"
  resumable_download "${BASE_URL}/md5.txt" "${dest}"
  local size; size="$(file_size_bytes "${dest}")"
  write_manifest_row "${MANIFEST}" "md5.txt" "${BASE_URL}/md5.txt" "${size}" "none" "none" "not_checked"
}

download_geneset() {
  local dest="${RAW_DIR}/GOPC.geneset.pep.fa.gz"
  local md5_file="${RAW_DIR}/md5.txt"
  if [ ! -f "${md5_file}" ]; then
    echo "md5.txt not present locally yet -- fetching it first to get the expected checksum." >&2
    download_md5
  fi
  local expected
  expected="$(grep -F 'GOPC.geneset.pep.fa.gz' "${md5_file}" | awk '{print $1}')"
  if [ -z "${expected}" ]; then
    echo "Could not find an md5 for GOPC.geneset.pep.fa.gz in ${md5_file} -- check its format hasn't changed upstream." >&2
    exit 2
  fi

  resumable_download "${BASE_URL}/GOPC.geneset.pep.fa.gz" "${dest}"
  local size; size="$(file_size_bytes "${dest}")"

  local status="verified_against_published"
  if ! verify_checksum md5 "${dest}" "${expected}"; then
    status="FAILED_verification"
  fi
  write_manifest_row "${MANIFEST}" "GOPC.geneset.pep.fa.gz" "${BASE_URL}/GOPC.geneset.pep.fa.gz" "${size}" "md5" "${expected}" "${status}"
  if [ "${status}" = "FAILED_verification" ]; then
    echo "Checksum verification FAILED -- see manifest and re-run to retry (curl -C - will resume, not restart)." >&2
    exit 1
  fi
}

case "${TARGET}" in
  md5) download_md5 ;;
  geneset) download_geneset ;;
  all) download_md5; download_geneset ;;
  *) echo "Unknown target '${TARGET}'. Use md5, geneset, or all." >&2; exit 2 ;;
esac
