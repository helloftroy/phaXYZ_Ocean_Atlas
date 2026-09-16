import csv

from phaatlas.pipeline import pathway_architecture as pa


def _write_metadata_tsv(path, family_id, rows):
    """rows: list of (target_id, genome, gtdb_genus, gtdb_species, study_id)."""
    fieldnames = [
        "target_id", "pha_family", "genome",
        "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species",
        "sample_id", "study_id", "latitude_degN", "longitude_degE",
        "ecosystem_type", "ecosystem_name", "ecosystem_compartment", "sample_source",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for target_id, genome, genus, species, study_id in rows:
            writer.writerow({
                "target_id": target_id, "pha_family": family_id, "genome": genome,
                "gtdb_domain": "Bacteria", "gtdb_phylum": "Pseudomonadota", "gtdb_class": "", "gtdb_order": "",
                "gtdb_family": "", "gtdb_genus": genus, "gtdb_species": species,
                "sample_id": f"{genome}_SAMPLE", "study_id": study_id,
                "latitude_degN": "1.0", "longitude_degE": "2.0",
                "ecosystem_type": "Marine", "ecosystem_name": "", "ecosystem_compartment": "", "sample_source": "",
            })


def test_short_code_derivation():
    assert pa.short_code("phaA") == "A"
    assert pa.short_code("phaJ") == "J"
    assert pa.short_code("phaR_regulator") == "RReg"
    assert pa.short_code("phaR_synthase") == "RSyn"


def test_load_genome_records_aggregates_across_families(tmp_path):
    _write_metadata_tsv(tmp_path / "phaA_unique_targets_with_metadata.tsv", "phaA", [
        ("T1", "G1", "Genus1", "Genus1 sp1", "STUDY1"),
        ("T2", "G1", "Genus1", "Genus1 sp1", "STUDY1"),  # second distinct phaA hit in same genome -> count 2
        ("T3", "G2", "Genus2", "Genus2 sp1", "STUDY2"),
    ])
    _write_metadata_tsv(tmp_path / "phaB_unique_targets_with_metadata.tsv", "phaB", [
        ("T4", "G1", "Genus1", "Genus1 sp1", "STUDY1"),  # G1 also has phaB -> architecture "AB"
        # G2 has no phaB -> architecture stays just "A"
    ])

    paths = pa.discover_metadata_files(tmp_path)
    assert len(paths) == 2

    genomes = pa.load_genome_records(paths)
    assert set(genomes) == {"G1", "G2"}
    assert genomes["G1"].family_counts() == {"phaA": 2, "phaB": 1}
    assert genomes["G2"].family_counts() == {"phaA": 1}
    assert genomes["G1"].taxonomy["gtdb_genus"] == "Genus1"


def test_architecture_label_uses_family_order():
    family_order = ["phaA", "phaB", "phaC", "phaJ"]
    assert pa.architecture_label({"phaC": 1, "phaA": 3}, family_order) == "AC"
    assert pa.architecture_label({"phaJ": 1, "phaB": 1, "phaA": 1}, family_order) == "ABJ"
    assert pa.architecture_label({}, family_order) == ""


def test_summarize_architectures_ranks_by_frequency(tmp_path):
    _write_metadata_tsv(tmp_path / "phaA_unique_targets_with_metadata.tsv", "phaA", [
        ("T1", "G1", "Genus1", "Genus1 sp1", "STUDY1"),
        ("T2", "G2", "Genus1", "Genus1 sp1", "STUDY1"),
        ("T3", "G3", "Genus2", "Genus2 sp1", "STUDY2"),
    ])
    _write_metadata_tsv(tmp_path / "phaB_unique_targets_with_metadata.tsv", "phaB", [
        ("T4", "G1", "Genus1", "Genus1 sp1", "STUDY1"),  # only G1 gets phaB -> architecture "AB"
    ])

    family_order = ["phaA", "phaB"]
    genomes = pa.load_genome_records(pa.discover_metadata_files(tmp_path))
    stats = pa.summarize_architectures(genomes, family_order, top_n=3)

    by_arch = {s.architecture: s for s in stats}
    assert by_arch["A"].n_genomes == 2  # G2, G3
    assert by_arch["AB"].n_genomes == 1  # G1
    # most common first
    assert stats[0].architecture == "A"
    assert stats[0].n_genomes == 2
    assert abs(by_arch["A"].pct_of_genomes - (2 / 3 * 100)) < 1e-6
    assert ("Genus1", 1) in by_arch["A"].top_genera or ("Genus2", 1) in by_arch["A"].top_genera


def test_write_genome_family_matrix_and_summary_round_trip(tmp_path):
    _write_metadata_tsv(tmp_path / "phaA_unique_targets_with_metadata.tsv", "phaA", [
        ("T1", "G1", "Genus1", "Genus1 sp1", "STUDY1"),
    ])
    family_order = ["phaA", "phaB"]
    genomes = pa.load_genome_records(pa.discover_metadata_files(tmp_path))

    matrix_path = tmp_path / "genome_family_matrix.tsv"
    n = pa.write_genome_family_matrix(genomes, family_order, matrix_path)
    assert n == 1
    with open(matrix_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["genome"] == "G1"
    assert rows[0]["n_phaA"] == "1"
    assert rows[0]["n_phaB"] == "0"
    assert rows[0]["architecture"] == "A"
    assert rows[0]["n_families_present"] == "1"

    summary_path = tmp_path / "architecture_summary.tsv"
    stats = pa.summarize_architectures(genomes, family_order)
    pa.write_architecture_summary(stats, summary_path)
    with open(summary_path, newline="") as f:
        summary_rows = list(csv.DictReader(f, delimiter="\t"))
    assert summary_rows[0]["architecture"] == "A"
    assert summary_rows[0]["n_genomes"] == "1"
