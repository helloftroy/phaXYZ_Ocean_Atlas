import csv

import httpx

from phaatlas.pipeline import ncbi_depth


def test_parse_depth_meters_common_formats():
    assert ncbi_depth.parse_depth_meters("0 m") == 0.0
    assert ncbi_depth.parse_depth_meters("200") == 200.0
    assert ncbi_depth.parse_depth_meters("surface") == 0.0
    assert ncbi_depth.parse_depth_meters("Surface water") == 0.0
    assert ncbi_depth.parse_depth_meters("10-20 m") == 15.0
    assert ncbi_depth.parse_depth_meters("10 to 20 m") == 15.0
    assert ncbi_depth.parse_depth_meters("-2 m") == -2.0


def test_parse_depth_meters_missing_values():
    assert ncbi_depth.parse_depth_meters(None) is None
    assert ncbi_depth.parse_depth_meters("") is None
    assert ncbi_depth.parse_depth_meters("not collected") is None
    assert ncbi_depth.parse_depth_meters("Not Applicable") is None
    assert ncbi_depth.parse_depth_meters("gibberish text with no number") is None


def test_depth_zone_bins():
    assert ncbi_depth.depth_zone(None) == ""
    assert ncbi_depth.depth_zone(-2.0) == "0-50m"  # clipped to 0
    assert ncbi_depth.depth_zone(0.0) == "0-50m"
    assert ncbi_depth.depth_zone(49.9) == "0-50m"
    assert ncbi_depth.depth_zone(50.0) == "50-200m"
    assert ncbi_depth.depth_zone(199.9) == "50-200m"
    assert ncbi_depth.depth_zone(200.0) == "200-1000m"
    assert ncbi_depth.depth_zone(999.9) == "200-1000m"
    assert ncbi_depth.depth_zone(1000.0) == "1000-4000m"
    assert ncbi_depth.depth_zone(3999.9) == "1000-4000m"
    assert ncbi_depth.depth_zone(4000.0) == ">4000m"
    assert ncbi_depth.depth_zone(11000.0) == ">4000m"


_SAMPLE_BIOSAMPLE_XML = """<?xml version="1.0" ?>
<BioSampleSet>
<BioSample access="public" id="1" accession="SAMN00000001">
  <Attributes>
    <Attribute attribute_name="depth" harmonized_name="depth" display_name="depth">0 m</Attribute>
    <Attribute attribute_name="lat_lon" harmonized_name="lat_lon">1.0 N 2.0 W</Attribute>
  </Attributes>
</BioSample>
<BioSample access="public" id="2" accession="SAMN00000002">
  <Attributes>
    <Attribute attribute_name="isolation_source" harmonized_name="isolation_source">sediment</Attribute>
  </Attributes>
</BioSample>
</BioSampleSet>"""


def test_fetch_sample_depths_batches_and_parses_xml():
    captured_params = []

    def handler(request: httpx.Request) -> httpx.Response:
        params = dict(httpx.QueryParams(request.url.query))
        captured_params.append(params)
        return httpx.Response(200, text=_SAMPLE_BIOSAMPLE_XML)

    client = httpx.Client(transport=httpx.MockTransport(handler))
    result = ncbi_depth.fetch_sample_depths(
        client, ["SAMN00000001", "SAMN00000002", "Unknown_biosample", ""], batch_size=10, sleep_seconds=0
    )

    assert result == {"SAMN00000001": "0 m", "SAMN00000002": None}
    assert len(captured_params) == 1  # both real accessions fit in one batch
    assert "SAMN00000001" in captured_params[0]["id"]
    assert "Unknown_biosample" not in captured_params[0]["id"]  # unresolvable marker filtered out


def test_fetch_sample_depths_respects_batch_size():
    calls = []

    def handler(request: httpx.Request) -> httpx.Response:
        calls.append(request.url.query)
        return httpx.Response(200, text='<?xml version="1.0"?><BioSampleSet></BioSampleSet>')

    client = httpx.Client(transport=httpx.MockTransport(handler))
    accessions = [f"SAMN{i:08d}" for i in range(5)]
    ncbi_depth.fetch_sample_depths(client, accessions, batch_size=2, sleep_seconds=0)
    assert len(calls) == 3  # 5 ids, batch_size=2 -> 2,2,1


def test_enrich_with_depth_end_to_end(tmp_path, monkeypatch):
    metadata_path = tmp_path / "phaC_unique_targets_with_metadata.tsv"
    with open(metadata_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["target_id", "genome", "biosample"])
        writer.writerow(["T1", "G1", "SAMN00000001"])
        writer.writerow(["T2", "G2", "SAMN00000002"])
        writer.writerow(["T3", "G3", "Unknown_biosample"])

    def fake_fetch(client, accessions, batch_size=100, sleep_seconds=0.35):
        return {"SAMN00000001": "0 m", "SAMN00000002": None}

    monkeypatch.setattr(ncbi_depth, "fetch_sample_depths", fake_fetch)
    monkeypatch.setattr(ncbi_depth, "_client", lambda: httpx.Client(transport=httpx.MockTransport(
        lambda request: httpx.Response(200, text='<?xml version="1.0"?><BioSampleSet></BioSampleSet>')
    )))

    out_path = tmp_path / "phaC_unique_targets_with_metadata_depth.tsv"
    stats = ncbi_depth.enrich_with_depth(metadata_path, out_path)

    assert stats["n_input_rows"] == 3
    assert stats["n_distinct_biosamples"] == 2
    assert stats["n_rows_with_parsed_depth"] == 1

    with open(out_path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    by_genome = {r["genome"]: r for r in rows}
    assert by_genome["G1"]["depth_raw"] == "0 m"
    assert by_genome["G1"]["depth_m"] == "0.0"
    assert by_genome["G1"]["depth_zone"] == "0-50m"
    assert by_genome["G2"]["depth_raw"] == ""
    assert by_genome["G2"]["depth_m"] == ""
    assert by_genome["G3"]["depth_zone"] == ""
