#!/usr/bin/env python3
"""Direct KOfam HMM scan for selected Thioglobus5 metabolic markers.

This avoids KOfamScan's Ruby/GNU-parallel wrapper while still using the
official KOfam profiles and ko_list score thresholds.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FASTA_DIR = ROOT / "thioglobus5_reannotation" / "proteins"
DEFAULT_WORK = ROOT / "thioglobus5_reannotation" / "kofam_marker_direct"
DEFAULT_SUMMARY = ROOT / "figures" / "thioglobus5_kofam_marker_summary.tsv"
DEFAULT_HITS = ROOT / "figures" / "thioglobus5_kofam_marker_hits.tsv"
DEFAULT_MISSING = ROOT / "figures" / "thioglobus5_kofam_marker_missing_profiles.tsv"

MARKERS = {
    "cbbM_rbcL": "K01601",
    "rbcS": "K01602",
    "prk": "K00855",
    "sqr": "K17218",
    "fccA": "K17230",
    "fccB": "K17229",
    "soxX": "K17223",
    "soxA": "K17222",
    "soxB": "K17224",
    "soxY": "K17226",
    "soxZ": "K17227",
    "soxC": "K17225",
    "soxD": "K22622",
    "dsrA": "K11180",
    "dsrB": "K11181",
    "dsrE": "K11182",
    "dsrF": "K11183",
    "dsrH": "K11184",
    "dsrC": "K11179",
    "aprA": "K00394",
    "aprB": "K00395",
    "sat": "K00958",
}
KO_TO_MARKER = {ko: marker for marker, ko in MARKERS.items()}
CBBCORE = {"cbbM_rbcL", "rbcS", "prk"}
SOX_CORE = {"soxX", "soxA", "soxB", "soxY", "soxZ"}


def read_ko_list(path: Path) -> dict[str, dict[str, str]]:
    meta: dict[str, dict[str, str]] = {}
    with path.open(errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header and header[0].lower() not in {"knum", "ko"}:
            handle.seek(0)
            header = [
                "knum",
                "threshold",
                "score_type",
                "profile_type",
                "f_measure",
                "nseq",
                "nseq_used",
                "alen",
                "mlen",
                "eff_nseq",
                "re_pos",
                "definition",
            ]
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            ko = row.get("knum") or row.get("ko")
            if ko in KO_TO_MARKER:
                meta[ko] = row
    return meta


def parse_domtbl(path: Path, ko_meta: dict[str, str]) -> list[dict[str, str]]:
    hits = []
    threshold = float(ko_meta.get("threshold") or 0)
    score_type = (ko_meta.get("score_type") or "full").lower()
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split(maxsplit=22)
            if len(parts) < 23:
                continue
            protein = parts[0]
            ko = parts[3]
            full_evalue = parts[6]
            full_score = float(parts[7])
            dom_evalue = parts[12]
            dom_score = float(parts[13])
            used_score = dom_score if score_type == "domain" else full_score
            if used_score < threshold:
                continue
            hits.append(
                {
                    "query": protein,
                    "ko": ko,
                    "marker": KO_TO_MARKER[ko],
                    "score_type": score_type,
                    "threshold": f"{threshold:g}",
                    "score": f"{used_score:g}",
                    "full_score": f"{full_score:g}",
                    "domain_score": f"{dom_score:g}",
                    "full_evalue": full_evalue,
                    "domain_evalue": dom_evalue,
                    "description": ko_meta.get("definition", ""),
                }
            )
    return hits


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        return {row["genome"]: row for row in csv.DictReader(handle, delimiter="\t")}


def find_profile(profiles_dir: Path, ko: str) -> Path | None:
    candidates = [
        profiles_dir / f"{ko}.hmm",
        profiles_dir / ko,
        profiles_dir / ko[1:3] / f"{ko}.hmm",
        profiles_dir / ko[1:3] / ko,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    for pattern in (f"{ko}.hmm", ko):
        matches = list(profiles_dir.rglob(pattern))
        if matches:
            return matches[0]
    return None


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kofam-db", required=True, type=Path)
    parser.add_argument("--hmmsearch-bin", default="hmmsearch")
    parser.add_argument("--fasta-dir", default=DEFAULT_FASTA_DIR, type=Path)
    parser.add_argument("--work-dir", default=DEFAULT_WORK, type=Path)
    parser.add_argument("--summary-out", default=DEFAULT_SUMMARY, type=Path)
    parser.add_argument("--hits-out", default=DEFAULT_HITS, type=Path)
    parser.add_argument("--missing-profiles-out", default=DEFAULT_MISSING, type=Path)
    parser.add_argument("--manifest", default=ROOT / "thioglobus5_reannotation" / "thioglobus5_manifest.tsv", type=Path)
    args = parser.parse_args()

    profiles = args.kofam_db / "profiles"
    ko_list = args.kofam_db / "ko_list"
    if not profiles.is_dir():
        raise SystemExit(f"Missing profiles directory: {profiles}")
    if not ko_list.is_file():
        raise SystemExit(f"Missing ko_list: {ko_list}")

    ko_meta = read_ko_list(ko_list)
    missing_meta = sorted(set(KO_TO_MARKER) - set(ko_meta))
    if missing_meta:
        raise SystemExit(f"Missing marker KOs in ko_list: {','.join(missing_meta)}")

    profile_by_ko: dict[str, Path] = {}
    missing_profiles = []
    for marker, ko in MARKERS.items():
        hmm = find_profile(profiles, ko)
        if hmm is None:
            missing_profiles.append(
                {
                    "marker": marker,
                    "ko": ko,
                    "profiles_dir": str(profiles),
                    "note": "KO is present in ko_list, but no KOfam HMM profile was found in this database release.",
                }
            )
        else:
            profile_by_ko[ko] = hmm

    args.missing_profiles_out.parent.mkdir(parents=True, exist_ok=True)
    with args.missing_profiles_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["marker", "ko", "profiles_dir", "note"])
        writer.writeheader()
        writer.writerows(missing_profiles)

    if missing_profiles:
        print(
            "Skipping marker KOs without KOfam profiles: "
            + ", ".join(f"{row['marker']}({row['ko']})" for row in missing_profiles)
        )

    manifest = read_manifest(args.manifest)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    all_hits: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for fasta in sorted(args.fasta_dir.glob("*.faa")):
        genome = fasta.stem
        genome_hits: list[dict[str, str]] = []
        genome_dir = args.work_dir / genome
        genome_dir.mkdir(parents=True, exist_ok=True)
        for marker, ko in MARKERS.items():
            hmm = profile_by_ko.get(ko)
            if hmm is None:
                continue
            domtbl = genome_dir / f"{ko}.domtbl"
            stdout = genome_dir / f"{ko}.hmmsearch.txt"
            cmd = [
                args.hmmsearch_bin,
                "--noali",
                "--domtblout",
                str(domtbl),
                str(hmm),
                str(fasta),
            ]
            with stdout.open("w") as out:
                subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, check=True)
            for hit in parse_domtbl(domtbl, ko_meta[ko]):
                hit["genome"] = genome
                hit["gtdb_species"] = manifest.get(genome, {}).get("gtdb_species", "")
                genome_hits.append(hit)
                all_hits.append(hit)

        present = {hit["marker"] for hit in genome_hits}
        row = {
            "genome": genome,
            "gtdb_species": manifest.get(genome, {}).get("gtdb_species", ""),
            "omdb_completeness_pct": manifest.get(genome, {}).get("omdb_completeness_pct", ""),
            "n_proteins": manifest.get(genome, {}).get("n_proteins", ""),
            "has_cbb_rbcLS_prk": yesno(CBBCORE <= present),
            "has_core_soxXABYZ": yesno(SOX_CORE <= present),
        }
        for marker in MARKERS:
            count = sum(1 for hit in genome_hits if hit["marker"] == marker)
            row[f"has_{marker}"] = yesno(count > 0)
            row[f"n_{marker}"] = str(count)
        summary_rows.append(row)

    summary_fields = [
        "genome",
        "gtdb_species",
        "omdb_completeness_pct",
        "n_proteins",
        "has_cbb_rbcLS_prk",
        "has_core_soxXABYZ",
    ]
    for marker in MARKERS:
        summary_fields += [f"has_{marker}", f"n_{marker}"]
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    with args.summary_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    hit_fields = [
        "genome",
        "gtdb_species",
        "marker",
        "ko",
        "query",
        "score_type",
        "threshold",
        "score",
        "full_score",
        "domain_score",
        "full_evalue",
        "domain_evalue",
        "description",
    ]
    with args.hits_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=hit_fields)
        writer.writeheader()
        writer.writerows(all_hits)

    print(f"Wrote {args.summary_out.relative_to(ROOT)}")
    print(f"Wrote {args.hits_out.relative_to(ROOT)}")
    print(f"Wrote {args.missing_profiles_out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
