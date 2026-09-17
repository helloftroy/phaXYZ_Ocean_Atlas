"""Deeper sequence clustering of a family's OMDB hits (default the coarse
"best_query" proxy in pipeline/depth_heatmap.py groups hits by whichever
NR95 reference protein they happened to match, not by similarity to each
other) -- this module supports actually clustering the DISCOVERED
sequences themselves at a chosen identity threshold (e.g. 50-70%), so
"does this protein cluster occupy a distinct environmental niche" is
asked about real sequence similarity groups, not reference-match groups.

Getting the sequences: rather than downloading OMDBv2.0_AA_G_NR100.faa.gz
again (~46GB, already avoided once for the cluster.tsv.gz join in
pipeline/omdb_metadata.py), this reuses the mmseqs target database
ALREADY BUILT for the search stage (build_gopc_target_db in
pipeline/gopc_search.py) via `mmseqs createsubdb --id-mode 1` (selects by
FASTA identifier, not internal DB key -- confirmed live against a small
local test target db) followed by `mmseqs convert2fasta`. Both mmseqs
calls and the identity-threshold clustering itself are cluster-side work
(need the actual target_db and a real mmseqs binary) -- see
cluster/run_omdb_extract_cluster_sequences.sbatch and
cluster/run_omdb_cluster_sequences.sbatch. This module's own Python code
only prepares the wanted-ID list beforehand and parses the resulting
cluster assignment TSV afterward -- both cheap, local, no mmseqs needed.
"""

from __future__ import annotations

import csv
from pathlib import Path


def collect_wanted_target_ids(metadata_paths: list[Path]) -> set[str]:
    """Distinct target_id values across one or more
    <family>_unique_targets_with_metadata[_depth].tsv files."""
    ids: set[str] = set()
    for path in metadata_paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                target_id = row.get("target_id", "")
                if target_id:
                    ids.add(target_id)
    return ids


def write_wanted_ids_file(target_ids: set[str], out_path: Path) -> int:
    """One target_id per line, sorted for a reproducible/diffable file --
    exactly the format `mmseqs createsubdb --id-mode 1` expects."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for target_id in sorted(target_ids):
            f.write(target_id + "\n")
    return len(target_ids)


def load_cluster_assignments(cluster_tsv_path: Path) -> dict[str, str]:
    """Parses mmseqs' own standard 2-column cluster TSV
    (representative \\t member, one row per member, from `mmseqs cluster`/
    `easy-cluster`'s own createtsv output) -- NOT the same 5-column
    CLUSTER/LENGTH/SIZE/REPRESENTATIVE/MEMBERS format OMDB's own
    OMDBv2.0_AA_G_NR100.cluster.tsv.gz uses (see
    pipeline/omdb_metadata.py's stream_target_cluster_genomes for that
    one) -- easy to mix these up since both are called "cluster.tsv" by
    convention. Returns {member_id: representative_id}, so every
    target_id (a member of exactly one cluster, including a
    singleton-cluster's own representative being its own sole member)
    maps to its cluster's representative id."""
    assignments: dict[str, str] = {}
    with open(cluster_tsv_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            representative, member = row[0], row[1]
            assignments[member] = representative
    return assignments
