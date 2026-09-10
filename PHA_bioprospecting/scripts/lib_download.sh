#!/usr/bin/env bash
# Shared helpers for download_gopc.sh / download_omdb.sh: resumable fetch
# (curl -C -), checksum verification, and a manifest row writer. Sourced,
# not run directly.
set -euo pipefail

# resumable_download <url> <dest_path>
# -C - resumes a partial download instead of restarting; -f fails loudly
# (rather than writing an HTML error page as if it were the file) on a
# non-2xx response; --retry handles transient network blips separately
# from a genuinely interrupted transfer, which -C - itself handles by
# re-running this same command.
resumable_download() {
  local url="$1" dest="$2"
  mkdir -p "$(dirname "${dest}")"
  echo "Downloading: ${url}"
  echo "        -->  ${dest}"
  curl -fL --retry 5 --retry-delay 5 -C - -o "${dest}" "${url}"
}

# compute_checksum <algo:md5|sha256> <file_path>  -> prints the hex digest
compute_checksum() {
  local algo="$1" file="$2"
  case "${algo}" in
    md5)
      if command -v md5sum >/dev/null 2>&1; then md5sum "${file}" | awk '{print $1}'
      else md5 -q "${file}"; fi  # macOS
      ;;
    sha256)
      if command -v sha256sum >/dev/null 2>&1; then sha256sum "${file}" | awk '{print $1}'
      else shasum -a 256 "${file}" | awk '{print $1}'; fi  # macOS
      ;;
    *) echo "compute_checksum: unknown algo '${algo}'" >&2; exit 2 ;;
  esac
}

# verify_checksum <algo> <file_path> <expected_hex>  -> exits 1 on mismatch
verify_checksum() {
  local algo="$1" file="$2" expected="$3"
  local actual
  actual="$(compute_checksum "${algo}" "${file}")"
  if [ "${actual}" != "${expected}" ]; then
    echo "CHECKSUM MISMATCH for ${file}" >&2
    echo "  expected (${algo}): ${expected}" >&2
    echo "  actual   (${algo}): ${actual}" >&2
    return 1
  fi
  echo "checksum OK (${algo}): ${file}"
  return 0
}

# write_manifest_row <manifest_tsv> <filename> <source_url> <size_bytes> <checksum_algo> <checksum> <checksum_status>
# checksum_status is one of: verified_against_published | computed_no_published_hash | not_checked
# Appends one row (writes the header first if the file is new). Never
# overwrites a prior row for the same filename -- appends a fresh one
# instead, so re-running a download after a source file changes upstream
# still leaves the old record for comparison rather than silently erasing it.
write_manifest_row() {
  local manifest="$1" filename="$2" source_url="$3" size_bytes="$4" checksum_algo="$5" checksum="$6" checksum_status="$7"
  mkdir -p "$(dirname "${manifest}")"
  if [ ! -f "${manifest}" ]; then
    printf 'filename\tsource_url\tsize_bytes\tchecksum_algo\tchecksum\tchecksum_status\tdownloaded_at\n' > "${manifest}"
  fi
  local downloaded_at
  downloaded_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${filename}" "${source_url}" "${size_bytes}" "${checksum_algo}" "${checksum}" "${checksum_status}" "${downloaded_at}" \
    >> "${manifest}"
  echo "manifest updated: ${manifest}"
}

file_size_bytes() {
  local file="$1"
  if command -v stat >/dev/null 2>&1; then
    stat -f%z "${file}" 2>/dev/null || stat -c%s "${file}" 2>/dev/null
  fi
}
