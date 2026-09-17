"""Prepare lightweight structure-prediction inputs for regional PhaC clusters.

Reads figures/phaC_regional_specialists_sequences.csv and writes:
  - one combined FASTA
  - one FASTA per cluster
  - metadata with cleaned sequence lengths and any residue cleanup notes

The exported protein sequences currently end in a single X, consistent with
a translated terminal stop/unknown marker. Structure predictors commonly
reject X, so this script trims only terminal X characters and records it.
"""
import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
IN_CSV = ROOT / "figures" / "phaC_regional_specialists_sequences.csv"
OUT = ROOT / "figures" / "structures" / "regional_specialists"
FASTA_DIR = OUT / "inputs"

AA = set("ACDEFGHIKLMNPQRSTVWY")


def slug(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return re.sub(r"_+", "_", value).strip("_")


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def clean_sequence(seq):
    raw = "".join(seq.split()).upper()
    trimmed = raw.rstrip("X")
    notes = []
    if len(trimmed) != len(raw):
        notes.append(f"trimmed_terminal_X:{len(raw) - len(trimmed)}")
    bad = sorted(set(trimmed) - AA)
    if bad:
        notes.append("noncanonical_after_cleanup:" + ",".join(bad))
    return raw, trimmed, ";".join(notes) if notes else "none"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(IN_CSV.open(newline="")))
    meta_path = OUT / "phaC_regional_specialists_structure_metadata.tsv"
    combined_path = OUT / "phaC_regional_specialists_for_structure.faa"

    with meta_path.open("w", newline="") as meta, combined_path.open("w") as combined:
        writer = csv.writer(meta, delimiter="\t")
        writer.writerow([
            "structure_id", "cluster_id", "dominant_genus", "n_genomes",
            "n_genomes_with_depth", "mean_depth_m", "latitude_degN",
            "longitude_degE", "raw_length", "clean_length", "cleanup_notes",
            "fasta_path",
        ])
        for row in rows:
            raw, clean, notes = clean_sequence(row["sequence"])
            short = row["cluster_id"][-12:]
            genus = slug(row["dominant_genus"])
            sid = f"{short}_{genus}"
            header = (
                f"{sid} cluster_id={row['cluster_id']} genus={row['dominant_genus']} "
                f"n_genomes={row['n_genomes']} mean_depth_m={row['mean_depth_m']} "
                f"lat={row['latitude_degN']} lon={row['longitude_degE']} cleanup={notes}"
            )
            fasta_path = FASTA_DIR / f"{sid}.faa"
            text = f">{header}\n{wrap(clean)}\n"
            fasta_path.write_text(text)
            combined.write(text)
            writer.writerow([
                sid, row["cluster_id"], row["dominant_genus"], row["n_genomes"],
                row["n_genomes_with_depth"], row["mean_depth_m"],
                row["latitude_degN"], row["longitude_degE"], len(raw), len(clean),
                notes, fasta_path.relative_to(ROOT),
            ])

    print(f"wrote {combined_path}")
    print(f"wrote {meta_path}")
    print(f"wrote {len(rows)} per-cluster FASTA files under {FASTA_DIR}")


if __name__ == "__main__":
    main()
