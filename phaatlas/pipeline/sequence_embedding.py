"""K-mer composition embedding of a family's extracted protein sequences
(see pipeline/sequence_clustering.py for how the FASTA is produced),
projected to 2D for exploratory sequence-space visualization -- "each
point is a protein", colorable downstream by sequence cluster, depth,
taxon, or pathway architecture.

PCA (scikit-learn, a core dependency) always works. UMAP is optional
(the `umap` extra, umap-learn -- heavier, pulls in numba) -- umap_2d
detects its absence and returns None rather than failing, so callers can
fall back to PCA-only.
"""

from __future__ import annotations

import csv
from itertools import product
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Returns [(id, sequence), ...] -- id is the first whitespace-delimited
    token of each header, matching mmseqs' own convention for how it
    assigns a sequence's "name" (and therefore its target_id)."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].split()[0] if line[1:].split() else line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def kmer_index(k: int, alphabet: str = AMINO_ACIDS) -> dict[str, int]:
    return {"".join(p): i for i, p in enumerate(product(alphabet, repeat=k))}


def kmer_composition_vector(sequence: str, k: int, index: dict[str, int]) -> np.ndarray:
    """Normalized k-mer frequency vector -- a k-mer containing a character
    outside `alphabet` (e.g. 'X' for an unresolved residue) is simply
    skipped rather than erroring, since OMDB sequences occasionally carry
    these."""
    vec = np.zeros(len(index), dtype=np.float64)
    n = 0
    for i in range(len(sequence) - k + 1):
        idx = index.get(sequence[i : i + k])
        if idx is not None:
            vec[idx] += 1
            n += 1
    if n > 0:
        vec /= n
    return vec


def build_feature_matrix(fasta_path: Path, k: int = 3) -> tuple[list[str], np.ndarray]:
    records = read_fasta(fasta_path)
    index = kmer_index(k)
    ids = [h for h, _ in records]
    if not records:
        return ids, np.zeros((0, len(index)))
    matrix = np.vstack([kmer_composition_vector(seq, k, index) for _, seq in records])
    return ids, matrix


def pca_2d(matrix: np.ndarray, random_state: int = 0) -> np.ndarray:
    """Always available (scikit-learn is a core dependency). Degenerates
    gracefully for very small inputs rather than raising: fewer than 2
    sequences, or fewer usable components than 2, pads with zero columns."""
    n = matrix.shape[0]
    if n < 2:
        return np.zeros((n, 2))
    n_components = min(2, n, matrix.shape[1])
    coords = PCA(n_components=n_components, random_state=random_state).fit_transform(matrix)
    if coords.shape[1] < 2:
        coords = np.pad(coords, ((0, 0), (0, 2 - coords.shape[1])))
    return coords


def umap_2d(matrix: np.ndarray, random_state: int = 0, n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray | None:
    """Returns None (not an exception) if umap-learn isn't installed, or
    if there are too few sequences for a meaningful embedding -- callers
    should fall back to pca_2d in either case."""
    try:
        import umap
    except ImportError:
        return None
    n = matrix.shape[0]
    if n < 4:
        return None
    reducer = umap.UMAP(n_components=2, random_state=random_state, n_neighbors=min(n_neighbors, n - 1), min_dist=min_dist)
    return reducer.fit_transform(matrix)


def collect_target_annotations(
    metadata_depth_path: Path,
    cluster_assignments: dict[str, str],
    genome_architecture: dict[str, str] | None = None,
) -> dict[str, dict]:
    """One annotation record per target_id, taken from the FIRST row seen
    for it (a target_id can appear against multiple genomes in the
    metadata file -- this picks a single representative genome's context
    for coloring a single embedding point, it does not attempt to
    reconcile disagreements across genomes)."""
    annotations: dict[str, dict] = {}
    with open(metadata_depth_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            target_id = row.get("target_id", "")
            if not target_id or target_id in annotations:
                continue
            genome = row.get("genome", "")
            annotations[target_id] = {
                "cluster_id": cluster_assignments.get(target_id, ""),
                "genome": genome,
                "genus": row.get("gtdb_genus", ""),
                "phylum": row.get("gtdb_phylum", ""),
                "study_id": row.get("study_id", ""),
                "depth_zone": row.get("depth_zone", ""),
                "architecture": (genome_architecture or {}).get(genome, ""),
            }
    return annotations


def write_embedding_tsv(
    ids: list[str],
    pca_coords: np.ndarray,
    umap_coords: np.ndarray | None,
    annotations: dict[str, dict],
    out_path: Path,
) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "target_id", "pca_x", "pca_y", "umap_x", "umap_y",
        "cluster_id", "genus", "phylum", "study_id", "depth_zone", "architecture",
    ]
    n = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for i, target_id in enumerate(ids):
            ann = annotations.get(target_id, {})
            row = {
                "target_id": target_id,
                "pca_x": f"{pca_coords[i][0]:.6f}", "pca_y": f"{pca_coords[i][1]:.6f}",
                "umap_x": f"{umap_coords[i][0]:.6f}" if umap_coords is not None else "",
                "umap_y": f"{umap_coords[i][1]:.6f}" if umap_coords is not None else "",
                "cluster_id": ann.get("cluster_id", ""),
                "genus": ann.get("genus", ""),
                "phylum": ann.get("phylum", ""),
                "study_id": ann.get("study_id", ""),
                "depth_zone": ann.get("depth_zone", ""),
                "architecture": ann.get("architecture", ""),
            }
            writer.writerow(row)
            n += 1
    return n
