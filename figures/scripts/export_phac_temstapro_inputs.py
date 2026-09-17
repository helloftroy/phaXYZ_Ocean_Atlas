"""Prepare all phaC target sequences for TemStaPro.

Input FASTA:
  data/temstapro_inputs/phaC_cluster_sequences.faa

Outputs:
  temstapro/phaC_temstapro_all.faa
  temstapro/chunks/phaC_temstapro_chunk_000.faa ...
  temstapro/phaC_temstapro_sequence_manifest.tsv

TemStaPro/ProtTrans can be sensitive to non-standard residues. We preserve
internal X characters, trim terminal X runs that look like translated stop or
unknown end markers, and replace rare noncanonical letters with X.
"""
import argparse
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_FASTA = ROOT / "data" / "temstapro_inputs" / "phaC_cluster_sequences.faa"
DEFAULT_OUT = ROOT / "temstapro"
AA = set("ACDEFGHIKLMNPQRSTVWYX")


def read_fasta(path):
    cur, desc, seq = None, "", []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    yield cur, desc, "".join(seq)
                desc = line[1:].strip()
                cur = desc.split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if cur:
            yield cur, desc, "".join(seq)


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def clean(seq):
    raw = "".join(seq.split()).upper()
    trimmed = raw.rstrip("X")
    notes = []
    if len(trimmed) != len(raw):
        notes.append(f"trimmed_terminal_X:{len(raw) - len(trimmed)}")
    fixed = []
    replaced = 0
    for aa in trimmed:
        if aa in AA:
            fixed.append(aa)
        else:
            fixed.append("X")
            replaced += 1
    if replaced:
        notes.append(f"noncanonical_to_X:{replaced}")
    return raw, "".join(fixed), ";".join(notes) if notes else "none"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-fasta", type=Path, default=DEFAULT_FASTA)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--chunk-size", type=int, default=1000)
    args = ap.parse_args()

    if not args.input_fasta.exists():
        raise SystemExit(
            f"Input FASTA not found: {args.input_fasta}\n"
            "This large FASTA is not stored in git. "
            "Copy or symlink phaC_cluster_sequences.faa into PHA_Ocean_Atlas, "
            "or pass its location explicitly, for example:\n"
            "  mkdir -p data/temstapro_inputs\n"
            "  ln -s /path/to/phaC_cluster_sequences.faa "
            "data/temstapro_inputs/phaC_cluster_sequences.faa\n"
            "  python3 figures/scripts/export_phac_temstapro_inputs.py "
            "--input-fasta /path/to/phaC_cluster_sequences.faa\n"
            "Default expected path:\n"
            "  PHA_Ocean_Atlas/data/temstapro_inputs/phaC_cluster_sequences.faa"
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    chunk_dir = args.out_dir / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    combined = args.out_dir / "phaC_temstapro_all.faa"
    manifest = args.out_dir / "phaC_temstapro_sequence_manifest.tsv"
    chunk_index = 0
    in_chunk = 0
    chunk_fh = None
    n = 0

    def open_chunk(idx):
        return (chunk_dir / f"phaC_temstapro_chunk_{idx:03d}.faa").open("w")

    with combined.open("w") as all_fh, manifest.open("w", newline="") as meta:
        writer = csv.writer(meta, delimiter="\t")
        writer.writerow([
            "target_id", "raw_length", "clean_length", "cleanup_notes",
            "chunk_index", "chunk_fasta", "original_header",
        ])
        for target_id, desc, seq in read_fasta(args.input_fasta):
            raw, cleaned, notes = clean(seq)
            if not cleaned:
                continue
            if chunk_fh is None or in_chunk >= args.chunk_size:
                if chunk_fh:
                    chunk_fh.close()
                chunk_fh = open_chunk(chunk_index)
                in_chunk = 0
                chunk_index += 1
            header = target_id
            record = f">{header}\n{wrap(cleaned)}\n"
            all_fh.write(record)
            chunk_fh.write(record)
            writer.writerow([
                target_id, len(raw), len(cleaned), notes, chunk_index - 1,
                f"chunks/phaC_temstapro_chunk_{chunk_index - 1:03d}.faa", desc,
            ])
            n += 1
            in_chunk += 1
    if chunk_fh:
        chunk_fh.close()
    print(f"wrote {n} sequences")
    print(f"combined FASTA: {combined}")
    print(f"chunks: {chunk_index}")
    print(f"manifest: {manifest}")


if __name__ == "__main__":
    main()
