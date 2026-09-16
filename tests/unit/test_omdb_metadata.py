import csv
import gzip
import json

import httpx
import pytest

from phaatlas.pipeline import omdb_metadata


def _write_cluster_tsv_gz(path, rows):
    """rows: list of (cluster, length, size, representative, members_list)."""
    with gzip.open(path, "wt", newline="") as f:
        f.write("CLUSTER\tLENGTH\tSIZE\tREPRESENTATIVE\tMEMBERS\n")
        for cluster, length, size, representative, members in rows:
            f.write(f"{cluster}\t{length}\t{size}\t{representative}\t{';'.join(members)}\n")


def test_parse_genome_from_member():
    assert omdb_metadata.parse_genome_from_member(
        "GARB21-1_SAMN12799101_MAG_00000001-scaffold_1_9"
    ) == "GARB21-1_SAMN12799101_MAG_00000001"


def test_parse_genome_from_member_rejects_unexpected_format():
    with pytest.raises(ValueError):
        omdb_metadata.parse_genome_from_member("not_a_real_member_id")


def test_stream_target_cluster_genomes_filters_and_dedupes(tmp_path):
    cluster_path = tmp_path / "cluster.tsv.gz"
    _write_cluster_tsv_gz(cluster_path, [
        # wanted: two members from the SAME genome (different scaffolds) -- must dedupe to 1 genome
        ("TGT1", 100, 2, "G1-scaffold_1_1", ["G1-scaffold_1_1", "G1-scaffold_2_5"]),
        # wanted: members spanning 3 distinct genomes, capped display at 2
        ("TGT2", 50, 3, "G1-scaffold_3_1", ["G1-scaffold_3_1", "G2-scaffold_1_1", "G3-scaffold_1_1"]),
        # NOT wanted -- must be skipped entirely (and must not raise, even with a differently-shaped id)
        ("TGT3", 10, 1, "whatever", ["whatever"]),
    ])

    result = omdb_metadata.stream_target_cluster_genomes(cluster_path, {"TGT1", "TGT2"}, max_genomes_per_target=2)

    assert set(result.keys()) == {"TGT1", "TGT2"}

    assert result["TGT1"].n_members == 2
    assert result["TGT1"].n_distinct_genomes == 1
    assert result["TGT1"].genomes == ["G1"]

    assert result["TGT2"].n_members == 3
    assert result["TGT2"].n_distinct_genomes == 3  # true count, uncapped
    assert result["TGT2"].genomes == ["G1", "G2"]  # capped at 2, sorted


def test_anchor_regex_batch_is_exact_not_substring():
    import re
    pattern = omdb_metadata._anchor_regex_batch(["abc", "de.f"])
    assert re.match(pattern, "abc")
    assert re.match(pattern, "de.f")
    assert not re.match(pattern, "ab")  # prefix must NOT match
    assert not re.match(pattern, "abcd")  # superstring must NOT match
    assert not re.match(pattern, "deXf")  # "." must be escaped, not treated as regex wildcard


def test_query_omdb_api_batches_and_sends_anchored_regex():
    captured_bodies = []

    def handler(request: httpx.Request) -> httpx.Response:
        body = json.loads(request.content)
        captured_bodies.append(body)
        # Echo back one row per id embedded in this batch's regex, so we can
        # confirm the response parsing path without a real server.
        value = body["search_value_pair"][0]["value"]
        ids_in_pattern = value.removeprefix("^(?:").removesuffix(")$").split("|")
        return httpx.Response(200, json={"data": [{"genome": i, "domain": "Bacteria"} for i in ids_in_pattern]})

    client = httpx.Client(transport=httpx.MockTransport(handler))
    ids = [f"G{i}" for i in range(5)]

    result = omdb_metadata._query_omdb_api(client, "genome-cols", "genome", ids, batch_size=2, sleep_seconds=0)

    assert len(captured_bodies) == 3  # 5 ids, batch_size=2 -> batches of 2,2,1
    assert captured_bodies[0]["search_value_pair"][0]["regex"] is True
    assert captured_bodies[0]["search_value_pair"][0]["value"].startswith("^(?:")
    assert captured_bodies[0]["search_value_pair"][0]["value"].endswith(")$")
    assert set(result.keys()) == set(ids)
    assert result["G0"]["domain"] == "Bacteria"


def test_enrich_unique_targets_end_to_end(tmp_path, monkeypatch):
    unique_targets_path = tmp_path / "phaX_unique_targets.tsv"
    with open(unique_targets_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "target_id", "pha_family", "best_query", "best_evalue", "best_bitscore",
            "best_pident", "best_qcov", "best_tcov", "n_reference_queries_matching",
        ])
        # TGT1: resolves to one genome, with full genome+sample metadata
        writer.writerow(["TGT1", "phaX", "Q1", "1e-10", "100", "80", "0.9", "0.9", "3"])
        # TGT2: resolves to a genome the fake API does NOT know about
        writer.writerow(["TGT2", "phaX", "Q2", "1e-8", "90", "70", "0.6", "0.6", "1"])
        # TGT3: not present in the cluster file at all
        writer.writerow(["TGT3", "phaX", "Q3", "1e-5", "50", "40", "0.3", "0.3", "1"])

    cluster_path = tmp_path / "cluster.tsv.gz"
    _write_cluster_tsv_gz(cluster_path, [
        ("TGT1", 100, 1, "G1-scaffold_1_1", ["G1-scaffold_1_1"]),
        ("TGT2", 100, 1, "G2-scaffold_1_1", ["G2-scaffold_1_1"]),
    ])

    def fake_fetch_genome_metadata(client, genome_ids, batch_size=100, api_base=None, sleep_seconds=0.3):
        assert set(genome_ids) == {"G1", "G2"}
        return {
            "G1": {
                "genome": "G1", "domain": "Bacteria", "phylum": "Pseudomonadota", "class_": "C",
                "order": "O", "family": "F", "genus": "Genus1", "species": "Genus1 sp1",
                "gtdb_version": "R226", "completeness": 95.0, "contamination": 1.0,
                "shushisample": "SAMPLE1", "sushistudy": "STUDY1",
            },
            # G2 deliberately absent -- simulates an unresolved genome
        }

    def fake_fetch_sample_metadata(client, sample_ids, batch_size=100, api_base=None, sleep_seconds=0.3):
        assert set(sample_ids) == {"SAMPLE1"}
        return {
            "SAMPLE1": {
                "sushisample": "SAMPLE1", "sushistudy": "STUDY1",
                "environment": "marine sediment metagenome", "ecosystem": "ocean",
                "ecosystem_type": "Marine", "ecosystem_name": "Marine benthic",
                "ecosystem_compartment": "Marine sediment", "sample_source": "Marine sediment",
                "latitude_degN": -76.55, "longitude_degE": -170.0,
                "biosample": "SAMN12799101", "bioproject": "PRJNA573088",
                "public_sample_link": "https://www.ncbi.nlm.nih.gov/biosample/SAMN12799101/",
            },
        }

    monkeypatch.setattr(omdb_metadata, "fetch_genome_metadata", fake_fetch_genome_metadata)
    monkeypatch.setattr(omdb_metadata, "fetch_sample_metadata", fake_fetch_sample_metadata)
    monkeypatch.setattr(omdb_metadata, "_client", lambda: httpx.Client(transport=httpx.MockTransport(
        lambda request: httpx.Response(200, json={"data": []})
    )))

    out_path = tmp_path / "phaX_unique_targets_with_metadata.tsv"
    summary = omdb_metadata.enrich_unique_targets(unique_targets_path, cluster_path, out_path)

    assert summary.n_input_rows == 3
    assert summary.n_target_ids == 3
    assert summary.n_targets_matched_in_cluster_file == 2  # TGT3 not in cluster file
    assert summary.n_distinct_genomes_needed == 2  # G1, G2
    assert summary.n_genomes_resolved == 1  # only G1
    assert summary.n_distinct_samples_needed == 1  # only from G1 (G2 unresolved has no sample)
    assert summary.n_samples_resolved == 1
    assert summary.n_output_rows == 3  # one row each for TGT1(->G1), TGT2(->G2), TGT3(no genome)

    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    by_target = {r["target_id"]: r for r in rows}

    tgt1 = by_target["TGT1"]
    assert tgt1["genome"] == "G1"
    assert tgt1["genome_metadata_found"] == "True"
    assert tgt1["gtdb_species"] == "Genus1 sp1"
    assert tgt1["sample_id"] == "SAMPLE1"
    assert tgt1["sample_metadata_found"] == "True"
    assert tgt1["latitude_degN"] == "-76.55"
    assert tgt1["ecosystem_compartment"] == "Marine sediment"
    # original unique_targets.tsv columns must survive untouched
    assert tgt1["best_query"] == "Q1"

    tgt2 = by_target["TGT2"]
    assert tgt2["genome"] == "G2"
    assert tgt2["genome_metadata_found"] == "False"
    assert tgt2["sample_metadata_found"] == "False"
    assert tgt2["gtdb_species"] == ""

    tgt3 = by_target["TGT3"]
    assert tgt3["genome"] == ""
    assert tgt3["genome_metadata_found"] == "False"
    assert tgt3["n_genomes_in_cluster_shown"] == "0"
