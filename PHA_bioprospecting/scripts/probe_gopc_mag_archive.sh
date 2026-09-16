#!/usr/bin/env bash
# Probe GOPC/GOMC MAG archive provenance without bulk extraction.
set -euo pipefail

ARCHIVE="${1:-PHA_bioprospecting/databases/GOPC/raw/43191.all_MAGs.tar.gz}"
HITS="${2:-hpc_results/phaB_unique_targets.tsv}"
OUT_DIR="${3:-hpc_results/gopc_mag_archive_probe}"

mkdir -p "${OUT_DIR}"

if [ ! -f "${ARCHIVE}" ]; then
  echo "Archive not found: ${ARCHIVE}" >&2
  echo "Place 43191.all_MAGs.tar.gz there, or pass its path as argument 1." >&2
  exit 2
fi

if [ ! -f "${HITS}" ]; then
  echo "GOPC hits table not found: ${HITS}" >&2
  exit 2
fi

FILELIST="${OUT_DIR}/43191_MAGs_filelist.txt"
IDS="${OUT_DIR}/representative_gopc_ids.tsv"
SUMMARY="${OUT_DIR}/archive_summary.txt"
MAPPING="${OUT_DIR}/gopc_mag_mapping_test.tsv"

tar -tzf "${ARCHIVE}" > "${FILELIST}.tmp"
head -100 "${FILELIST}.tmp" > "${OUT_DIR}/43191_MAGs_filelist_head100.txt"
mv "${FILELIST}.tmp" "${FILELIST}"

awk 'NR == 1 { next }
     $1 ~ /^[0-9]+_MG_FD_gene_id_[0-9]+$/ {
       split($1, a, "_gene_id_");
       if (++seen_prefix[a[1]] <= 3 && total < 20) {
         print $1 "\t" a[1];
         total++;
       }
     }' "${HITS}" > "${IDS}"

{
  echo "archive	${ARCHIVE}"
  echo "file_count	$(wc -l < "${FILELIST}" | tr -d ' ')"
  echo
  echo "[top_level_paths]"
  awk -F/ '{ print $1 }' "${FILELIST}" | sort | uniq -c | sort -nr | head -50
  echo
  echo "[extensions]"
  awk '
    {
      name=$0
      sub(/^.*\//, "", name)
      ext="(none)"
      if (name ~ /\.[^.]+$/) {
        ext=name
        sub(/^.*\./, ".", ext)
      }
      print ext
    }' "${FILELIST}" | sort | uniq -c | sort -nr
  echo
  echo "[filelist_mg_fd_examples]"
  grep -E '_MG_FD' "${FILELIST}" | head -50 || true
  echo
  echo "[candidate_metadata_files]"
  grep -Ei 'readme|metadata|meta|tax|taxonomy|mapping|map|sample|biosample|sra|run|\.tsv$|\.csv$|\.gff($|\.gz$)|\.gbk($|\.gz$)|\.faa($|\.gz$)|\.fna($|\.gz$)' "${FILELIST}" | head -200 || true
} > "${SUMMARY}"

{
  printf 'gopc_gene_id\tgopc_prefix\tfound_in_mag_archive\tmag_id\tcontig_id\tprotein_id_in_mag\tsample_id\tsra_run\ttaxonomy\tmetadata_source\tnotes\n'
  while IFS=$'\t' read -r gene prefix; do
    prefix_hits="$(grep -F "${prefix}" "${FILELIST}" | head -20 || true)"
    gene_hits="$(grep -F "${gene}" "${FILELIST}" | head -20 || true)"
    if [ -n "${gene_hits}" ] || [ -n "${prefix_hits}" ]; then
      first_path="$(printf '%s\n%s\n' "${gene_hits}" "${prefix_hits}" | sed '/^$/d' | head -1)"
      printf '%s\t%s\tfilelist_match\t\t\t\t\t\t\t43191_MAGs_filelist.txt\tmatched path: %s\n' "${gene}" "${prefix}" "${first_path}"
    else
      printf '%s\t%s\tno\t\t\t\t\t\t\t43191_MAGs_filelist.txt\tnot found among MAG archive paths; this does not rule out unbinned GOPC contigs\n' "${gene}" "${prefix}"
    fi
  done < "${IDS}"
} > "${MAPPING}"

while IFS=$'\t' read -r gene prefix; do
  grep -F "${prefix}" "${FILELIST}" | head -20 > "${OUT_DIR}/filelist_matches_${prefix}.txt" || true
  grep -F "${gene}" "${FILELIST}" | head -20 > "${OUT_DIR}/filelist_matches_${gene}.txt" || true
done < "${IDS}"

echo "Wrote:"
echo "  ${FILELIST}"
echo "  ${SUMMARY}"
echo "  ${IDS}"
echo "  ${MAPPING}"
echo
echo "Next manual step, if a matched FAA/GFF/FNA path appears:"
echo "  tar -xOf ${ARCHIVE} path/to/matched/file | head -100"
