"""Exact-dedup + per-family MMseqs2 clustering at 95% identity / 90%
coverage, writing cluster95_id/cluster95_representative/cluster95_size back
onto family_assignment. Never deletes anything -- only ever ADD COLUMN
migrations (db/session.py) and in-place UPDATEs of those three columns.

Two-stage reduction, and why: MMseqs2 clusters the exact-duplicate-collapsed
set (one entry per unique sequence_sha256 within a family), not every row,
so a family with heavy exact-sequence redundancy (the same PhaC sequence
deposited under five accessions) doesn't waste alignment time re-comparing
identical sequences to each other. Every actual protein/family row sharing
that hash is then stamped with the SAME cluster95_id afterward, so
cluster95_size reflects the real number of database entries in the cluster,
not just the number of distinct sequences MMseqs2 actually saw.
"""

from __future__ import annotations

import datetime as dt
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from phaatlas.db.models import FamilyAssignment, Protein, RetrievalRun

EVIDENCE_TIER_RANK = {
    "GOLD_EXPERIMENTAL": 4,
    "CURATED_REVIEWED": 3,
    "ANNOTATED_UNREVIEWED": 2,
    "CANDIDATE_AMBIGUOUS": 1,
}


class MMseqsNotFoundError(RuntimeError):
    pass


@dataclass
class FamilyMember:
    protein_id: str
    sequence: str
    sequence_sha256: str
    evidence_tier: str


def _fetch_family_members(session: Session, family_id: str) -> list[FamilyMember]:
    rows = session.execute(
        select(Protein.protein_id, Protein.sequence, Protein.sequence_sha256, FamilyAssignment.evidence_tier)
        .join(FamilyAssignment, FamilyAssignment.protein_id == Protein.protein_id)
        .where(FamilyAssignment.family_id == family_id)
        .where(Protein.sequence.is_not(None))
        .where(Protein.sequence != "")
    ).all()
    return [FamilyMember(protein_id=r[0], sequence=r[1], sequence_sha256=r[2], evidence_tier=r[3]) for r in rows]


def _pick_canonical_per_hash(
    members: list[FamilyMember],
) -> tuple[dict[str, list[str]], dict[str, str], dict[str, str]]:
    """Groups members by sequence_sha256. Returns:
    - hash_to_protein_ids: every protein_id sharing that exact sequence
    - hash_to_canonical_seq: the one sequence text to feed MMseqs2 per hash
      (all members of a hash group are byte-identical by definition, so any
      one of them works -- picked deterministically: best evidence_tier,
      then lowest protein_id, so re-running is reproducible).
    """
    by_hash: dict[str, list[FamilyMember]] = {}
    for m in members:
        by_hash.setdefault(m.sequence_sha256, []).append(m)

    hash_to_protein_ids: dict[str, list[str]] = {}
    hash_to_canonical: dict[str, FamilyMember] = {}
    for h, group in by_hash.items():
        group_sorted = sorted(group, key=lambda m: (-EVIDENCE_TIER_RANK.get(m.evidence_tier, 0), m.protein_id))
        hash_to_canonical[h] = group_sorted[0]
        hash_to_protein_ids[h] = [m.protein_id for m in group]

    hash_to_canonical_id = {h: m.protein_id for h, m in hash_to_canonical.items()}
    hash_to_canonical_seq = {h: m.sequence for h, m in hash_to_canonical.items()}
    return hash_to_protein_ids, hash_to_canonical_id, hash_to_canonical_seq


def _write_fasta(canonical_seqs: dict[str, str], canonical_ids: dict[str, str], path: Path) -> None:
    with open(path, "w") as f:
        for h, seq in canonical_seqs.items():
            protein_id = canonical_ids[h]
            f.write(f">{protein_id}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def _mmseqs_version(mmseqs_bin: str) -> str:
    try:
        out = subprocess.run([mmseqs_bin, "version"], capture_output=True, text=True, check=True)
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def run_mmseqs_easy_cluster(
    fasta_path: Path,
    work_dir: Path,
    min_seq_id: float,
    min_cov: float,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
) -> dict[str, str]:
    """Runs `mmseqs easy-cluster` and returns {member_protein_id: representative_protein_id}
    over the DEDUPED (canonical-per-hash) input set only -- caller expands
    to every actual protein_id sharing that hash afterward.
    """
    if shutil.which(mmseqs_bin) is None:
        raise MMseqsNotFoundError(
            f"'{mmseqs_bin}' not found on PATH. Install MMseqs2 (conda: `conda install -c bioconda -c conda-forge "
            f"mmseqs2`, or a static binary from https://github.com/soedinglab/MMseqs2/releases) and re-run, "
            f"or pass --mmseqs-bin with the full path."
        )
    out_prefix = work_dir / "result"
    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        mmseqs_bin, "easy-cluster",
        str(fasta_path), str(out_prefix), str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(min_cov),
        "--cov-mode", "0",
    ]
    if threads:
        cmd += ["--threads", str(threads)]

    subprocess.run(cmd, check=True, capture_output=True, text=True)

    cluster_tsv = Path(f"{out_prefix}_cluster.tsv")
    mapping: dict[str, str] = {}
    with open(cluster_tsv) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rep, member = line.split("\t")
            mapping[member] = rep
    return mapping


def cluster_family(
    session: Session,
    family_id: str,
    min_seq_id: float = 0.95,
    min_cov: float = 0.90,
    mmseqs_bin: str = "mmseqs",
    threads: int | None = None,
) -> dict:
    """Clusters one family end-to-end and writes cluster95_* back onto every
    family_assignment row for it that has a sequence. Returns a summary dict.
    """
    members = _fetch_family_members(session, family_id)
    summary = {"family_id": family_id, "members_with_sequence": len(members), "unique_sequences": 0, "clusters": 0}
    if not members:
        return summary

    hash_to_protein_ids, hash_to_canonical_id, hash_to_canonical_seq = _pick_canonical_per_hash(members)
    summary["unique_sequences"] = len(hash_to_canonical_id)

    with tempfile.TemporaryDirectory(prefix=f"cluster95_{family_id}_") as tmp:
        work_dir = Path(tmp)
        fasta_path = work_dir / "dedup.faa"
        _write_fasta(hash_to_canonical_seq, hash_to_canonical_id, fasta_path)

        member_to_rep = run_mmseqs_easy_cluster(
            fasta_path, work_dir, min_seq_id=min_seq_id, min_cov=min_cov, mmseqs_bin=mmseqs_bin, threads=threads
        )
        version = _mmseqs_version(mmseqs_bin)

    # member_to_rep is keyed by the CANONICAL protein_id we fed MMseqs2
    # (one per hash) -- expand back to every real protein_id sharing that
    # hash, and count cluster size across that full expansion.
    rep_counts: dict[str, int] = {}
    protein_id_to_rep: dict[str, str] = {}
    for h, canonical_id in hash_to_canonical_id.items():
        rep = member_to_rep.get(canonical_id)
        if rep is None:
            # Shouldn't happen (every input sequence appears in its own
            # cluster's tsv at minimum, even as a singleton) -- skip rather
            # than crash the whole family if MMseqs2's output is ever short.
            continue
        for protein_id in hash_to_protein_ids[h]:
            protein_id_to_rep[protein_id] = rep

    for rep in protein_id_to_rep.values():
        rep_counts[rep] = rep_counts.get(rep, 0) + 1
    summary["clusters"] = len(rep_counts)

    assignments = (
        session.query(FamilyAssignment)
        .filter(FamilyAssignment.family_id == family_id, FamilyAssignment.protein_id.in_(protein_id_to_rep.keys()))
        .all()
    )
    now = dt.datetime.now(dt.timezone.utc)
    for fa in assignments:
        rep = protein_id_to_rep[fa.protein_id]
        fa.cluster95_id = f"{family_id}:{rep}"
        fa.cluster95_representative = rep
        fa.cluster95_size = rep_counts[rep]
        fa.updated_at = now

    run = RetrievalRun(
        database="MMseqs2",
        database_release=version,
        family_id=family_id,
        query=f"easy-cluster --min-seq-id {min_seq_id} -c {min_cov} --cov-mode 0",
        number_returned=summary["clusters"],
        software_version="phaatlas/0.1",
        status="completed",
        finished_at=now,
    )
    session.add(run)
    session.commit()
    return summary
