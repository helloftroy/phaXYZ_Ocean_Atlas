"""Resolve phaC target IDs to all OMDB genomes via the NR100 cluster file.

This builds an uncapped numerator for "% of all genomes with phaC".
The ordinary *_with_metadata.tsv files may cap how many genomes are shown
per NR100 target cluster; this script streams the full cluster.tsv.gz once
and emits every distinct genome represented in each retained phaC NR100
cluster.

Important: the original MMseqs search was run with permissive coverage
settings so fragmented hits were not lost at search time. For prevalence,
do not expand every weak target. By default this script keeps only target
clusters whose best phaC alignment has evalue <= 1e-20 and both query and
target coverage >= 0.5.

Default output:
  data/all_genomes/phaC_all_genomes_from_nr100_clusters.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_PHAC_TARGETS = ROOT / "PHA_bioprospecting" / "omdb_search" / "results" / "phaC_unique_targets.tsv"
DEFAULT_CLUSTER_TSV = ROOT / "PHA_bioprospecting" / "databases" / "OMDBv2" / "OMDBv2.0_AA_G_NR100.cluster.tsv.gz"
DEFAULT_OUT = ROOT / "data" / "all_genomes" / "phaC_all_genomes_from_nr100_clusters.tsv"


def parse_genome_from_member(member_id: str) -> str:
    marker = "-scaffold_"
    idx = member_id.find(marker)
    if idx == -1:
        return ""
    return member_id[:idx]


def get_float(row: dict[str, str], col: str, default: float) -> float:
    raw = row.get(col, "")
    if raw == "":
        return default
    return float(raw)


def load_target_ids(
    path: Path,
    max_evalue: float,
    min_qcov: float,
    min_tcov: float,
    min_pident: float,
    min_bitscore: float,
) -> tuple[set[str], int]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "target_id" not in reader.fieldnames:
            raise SystemExit(f"{path} does not have a target_id column")
        out = set()
        n = 0
        for row in reader:
            if not row.get("target_id"):
                continue
            n += 1
            if (
                get_float(row, "best_evalue", 1.0) <= max_evalue
                and get_float(row, "best_qcov", 0.0) >= min_qcov
                and get_float(row, "best_tcov", 0.0) >= min_tcov
                and get_float(row, "best_pident", 0.0) >= min_pident
                and get_float(row, "best_bitscore", 0.0) >= min_bitscore
            ):
                out.add(row["target_id"])
        return out, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phac-targets", type=Path, default=DEFAULT_PHAC_TARGETS)
    ap.add_argument("--cluster-tsv", type=Path, default=DEFAULT_CLUSTER_TSV)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--max-evalue", type=float, default=1e-20)
    ap.add_argument("--min-qcov", type=float, default=0.5)
    ap.add_argument("--min-tcov", type=float, default=0.5)
    ap.add_argument("--min-pident", type=float, default=0.0)
    ap.add_argument("--min-bitscore", type=float, default=0.0)
    args = ap.parse_args()

    if not args.phac_targets.exists():
        raise SystemExit(f"phaC target table not found: {args.phac_targets}")
    if not args.cluster_tsv.exists():
        raise SystemExit(
            f"OMDB NR100 cluster file not found: {args.cluster_tsv}\n"
            "Download it on the cluster with:\n"
            "  cd PHA_bioprospecting/scripts\n"
            "  ./download_omdb.sh nr100-clusters"
        )

    wanted, n_input_targets = load_target_ids(
        args.phac_targets,
        max_evalue=args.max_evalue,
        min_qcov=args.min_qcov,
        min_tcov=args.min_tcov,
        min_pident=args.min_pident,
        min_bitscore=args.min_bitscore,
    )
    print(f"input phaC target IDs: {n_input_targets:,}")
    print(
        "retained phaC target IDs after filters "
        f"(e<={args.max_evalue:g}, qcov>={args.min_qcov:g}, "
        f"tcov>={args.min_tcov:g}, pident>={args.min_pident:g}, "
        f"bits>={args.min_bitscore:g}): {len(wanted):,}"
    )
    if not wanted:
        raise SystemExit("No target IDs passed filters.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    matched_targets = 0
    output_pairs = set()
    with gzip.open(args.cluster_tsv, "rt", newline="") as f, args.out.open("w", newline="") as out_f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            cluster_idx = header.index("CLUSTER")
            members_idx = header.index("MEMBERS")
        except ValueError:
            cluster_idx = 0
            members_idx = 4
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow([
            "target_id", "genome", "max_evalue", "min_qcov", "min_tcov",
            "min_pident", "min_bitscore",
        ])
        for line in f:
            first_tab = line.find("\t")
            cluster_id = line[:first_tab]
            if cluster_id not in wanted:
                continue
            fields = line.rstrip("\n").split("\t")
            target_id = fields[cluster_idx]
            genomes = {parse_genome_from_member(m) for m in fields[members_idx].split(";")}
            genomes.discard("")
            for genome in sorted(genomes):
                pair = (target_id, genome)
                if pair not in output_pairs:
                    writer.writerow([
                        target_id, genome, args.max_evalue, args.min_qcov,
                        args.min_tcov, args.min_pident, args.min_bitscore,
                    ])
                    output_pairs.add(pair)
            matched_targets += 1

    distinct_genomes = len({genome for _, genome in output_pairs})
    print(f"matched target IDs in cluster file: {matched_targets:,}")
    print(f"target-genome rows: {len(output_pairs):,}")
    print(f"distinct phaC genomes: {distinct_genomes:,}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
