"""Merge TemStaPro chunk outputs and join phaC metadata + WOA temperature."""
import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_BASE = ROOT / "temstapro"
DEFAULT_INPUT_DIR = ROOT / "data" / "temstapro_inputs"
DEFAULT_MANIFEST = DEFAULT_BASE / "phaC_temstapro_sequence_manifest.tsv"
DEFAULT_WOA = DEFAULT_INPUT_DIR / "phaC_genomes_woa23_annual_temperature.tsv"
DEFAULT_META = DEFAULT_INPUT_DIR / "phaC_unique_targets_with_metadata_depth.tsv"
DEFAULT_OUT = DEFAULT_BASE / "phaC_temstapro_predictions_with_metadata.tsv"


def load_by_key(path, key, delimiter="\t"):
    out = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f, delimiter=delimiter):
            out[row[key]] = row
    return out


def require_file(path, label):
    if not path.exists():
        raise SystemExit(
            f"{label} not found: {path}\n"
            "Expected TemStaPro inputs under PHA_Ocean_Atlas/data/temstapro_inputs "
            "and generated outputs under PHA_Ocean_Atlas/temstapro. "
            "Use the corresponding --*-path option if this file lives elsewhere."
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--woa", type=Path, default=DEFAULT_WOA)
    ap.add_argument("--metadata", type=Path, default=DEFAULT_META)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    require_file(args.manifest, "Sequence manifest")
    require_file(args.woa, "WOA temperature table")
    require_file(args.metadata, "phaC metadata table")

    manifest = load_by_key(args.manifest, "target_id")
    genome_woa = load_by_key(args.woa, "genome")

    target_meta = {}
    with args.metadata.open(newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            target_meta.setdefault(row["target_id"], row)

    pred_files = sorted((args.base_dir / "temstapro_outputs").glob("phaC_temstapro_chunk_*.tsv"))
    if not pred_files:
        raise SystemExit(f"No TemStaPro outputs found under {args.base_dir / 'temstapro_outputs'}")

    pred_rows = []
    pred_fields = None
    for p in pred_files:
        with p.open(newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            if pred_fields is None:
                pred_fields = r.fieldnames
            for row in r:
                pred_rows.append(row)

    # TemStaPro typically uses the FASTA header as the first column; be robust
    # to minor header-name changes.
    id_col = None
    for c in pred_fields:
        if c.lower() in {"protein", "protein_id", "seq", "sequence", "fasta_header", "id", "header"}:
            id_col = c
            break
    if id_col is None:
        id_col = pred_fields[0]

    meta_cols = [
        "genome", "sample_id", "study_id", "biosample", "bioproject",
        "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
        "gtdb_family", "gtdb_genus", "gtdb_species", "latitude_degN",
        "longitude_degE", "depth_m", "depth_zone",
    ]
    woa_cols = [
        "woa23_temp_annual_degC", "woa23_lat", "woa23_lon",
        "woa23_horizontal_distance_km", "woa23_depth_method",
        "woa23_match_status",
    ]
    out_fields = ["target_id"] + pred_fields + [
        "raw_length", "clean_length", "cleanup_notes",
    ] + meta_cols + woa_cols

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=out_fields, extrasaction="ignore")
        w.writeheader()
        for pred in pred_rows:
            tid = pred[id_col].split()[0]
            meta = target_meta.get(tid, {})
            woa = genome_woa.get(meta.get("genome", ""), {})
            man = manifest.get(tid, {})
            row = {"target_id": tid}
            row.update(pred)
            row.update({k: man.get(k, "") for k in ["raw_length", "clean_length", "cleanup_notes"]})
            row.update({k: meta.get(k, "") for k in meta_cols})
            row.update({k: woa.get(k, "") for k in woa_cols})
            w.writerow(row)
    print(f"merged {len(pred_rows)} predictions -> {args.out}")


if __name__ == "__main__":
    main()
