#!/usr/bin/env bash
# Downloads OMDBv2 (Ocean Microbiomics Database v2, Sunagawa Lab / ETH
# Zurich -- same long-running project as Paoli, Ruscheweyh, Forneris et
# al., "Biosynthetic potential of the global ocean microbiome", Nature
# 607:111-118 (2022), doi:10.1038/s41586-022-04862-3; v2.0 is a later data
# release of the same portal, not a separately published paper as far as
# could be confirmed). Files are served directly from ETH's own file
# server (no login), confirmed to support HTTP Range requests (resumable).
#
# Usage:
#   ./download_omdb.sh data              # OMDBv2.0_data.tsv.gz -- ~7.5MB, genome/sample/study links. Safe anywhere.
#   ./download_omdb.sh nr100             # OMDBv2.0_AA_G_NR100.faa.gz -- ~49GB (100%-identity dedup AA catalog). Cluster only.
#   ./download_omdb.sh nr100-clusters    # OMDBv2.0_AA_G_NR100.cluster.tsv.gz -- ~4.2GB (cluster->genome membership). Cluster only.
#   ./download_omdb.sh all               # all three
#
# Only OMDBv2.0_data.tsv.gz has a published checksum (MD5, from the portal's
# own suppl_info page). The two large NR100 files have none published --
# this script computes and records our OWN sha256 for them instead, so at
# least future re-downloads/transfers of OUR copy can be checked for
# corruption, even without an upstream hash to compare against.
set -euo pipefail
cd "$(dirname "$0")"
source ./lib_download.sh

BASE_URL="https://sunagawalab.ethz.ch/share/microbiomics/ocean/db/2.0/data/catalogs"
OUT_DIR="../databases/OMDBv2"
MANIFEST="../databases/OMDBv2/download_manifest.tsv"

TARGET="${1:-}"
if [ -z "${TARGET}" ]; then
  echo "Usage: $0 {data|nr100|nr100-clusters|all}" >&2
  exit 2
fi

download_data() {
  local dest="${OUT_DIR}/OMDBv2.0_data.tsv.gz"
  local url="${BASE_URL}/OMDBv2.0_data.tsv.gz"
  local expected="c1b5f14c9b7899f7300ccf41e62f8681"  # published MD5, from https://omdb.microbiomics.io/repository/ocean/suppl_info

  resumable_download "${url}" "${dest}"
  local size; size="$(file_size_bytes "${dest}")"

  local status="verified_against_published"
  if ! verify_checksum md5 "${dest}" "${expected}"; then
    status="FAILED_verification"
  fi
  write_manifest_row "${MANIFEST}" "OMDBv2.0_data.tsv.gz" "${url}" "${size}" "md5" "${expected}" "${status}"
  if [ "${status}" = "FAILED_verification" ]; then
    echo "Checksum verification FAILED -- see manifest and re-run to retry." >&2
    exit 1
  fi
}

download_no_published_hash() {
  local filename="$1"
  local dest="${OUT_DIR}/${filename}"
  local url="${BASE_URL}/OMDBv2.0_AA_G_NR100/${filename}"

  resumable_download "${url}" "${dest}"
  local size; size="$(file_size_bytes "${dest}")"

  echo "No published checksum for ${filename} -- computing sha256 of our own copy for future reference..."
  local our_hash; our_hash="$(compute_checksum sha256 "${dest}")"
  echo "sha256: ${our_hash}"
  write_manifest_row "${MANIFEST}" "${filename}" "${url}" "${size}" "sha256" "${our_hash}" "computed_no_published_hash"
}

case "${TARGET}" in
  data) download_data ;;
  nr100) download_no_published_hash "OMDBv2.0_AA_G_NR100.faa.gz" ;;
  nr100-clusters) download_no_published_hash "OMDBv2.0_AA_G_NR100.cluster.tsv.gz" ;;
  all)
    download_data
    download_no_published_hash "OMDBv2.0_AA_G_NR100.faa.gz"
    download_no_published_hash "OMDBv2.0_AA_G_NR100.cluster.tsv.gz"
    ;;
  *) echo "Unknown target '${TARGET}'. Use data, nr100, nr100-clusters, or all." >&2; exit 2 ;;
esac
