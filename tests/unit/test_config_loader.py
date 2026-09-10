from phaatlas.config_loader import (
    DEFAULT_FAMILY_CONFIG_PATH,
    load_family_definitions,
    load_monomer_vocab,
    load_polymer_vocab,
)

EXPECTED_FAMILY_IDS = {
    "phaA",
    "phaB",
    "phaC",
    "phaE",
    "phaR_synthase",
    "phaJ",
    "phaG",
    "phaZ",
    "phaY",
    "phaP",
    "phaF",
    "phaI",
    "phaM",
    "phaD",
    "phaR_regulator",
    "phaQ",
}


def test_loads_all_16_families_from_spec():
    families = load_family_definitions()
    ids = {f.family_id for f in families}
    assert ids == EXPECTED_FAMILY_IDS


def test_phaR_synthase_and_regulator_are_distinct_families_sharing_gene_phaR():
    families = {f.family_id: f for f in load_family_definitions()}
    synth = families["phaR_synthase"]
    reg = families["phaR_regulator"]
    assert synth is not reg
    assert "phaR" in synth.genes
    assert "phaR" in reg.genes
    # each needs its own disambiguation terms to ever be told apart
    assert synth.disambiguation_terms
    assert reg.disambiguation_terms
    assert set(synth.disambiguation_terms).isdisjoint(reg.disambiguation_terms)


def test_only_phaC_allows_broad_uniprot_search():
    families = load_family_definitions()
    broad = [f.family_id for f in families if f.broad_uniprot_search]
    assert broad == ["phaC"]


def test_known_ec_mappings_from_spec():
    families = {f.family_id: f for f in load_family_definitions()}
    assert families["phaC"].ec_brenda == "2.3.1.304"
    assert families["phaA"].ec_brenda == "2.3.1.9"
    assert families["phaB"].ec_brenda == "1.1.1.36"
    assert families["phaJ"].ec_brenda == "4.2.1.119"


def test_ambiguous_ec_families_require_pha_context():
    families = {f.family_id: f for f in load_family_definitions()}
    assert families["phaA"].require_pha_context is True
    assert families["phaJ"].require_pha_context is True
    # phaC's EC is clean/specific enough to not need the extra gate
    assert families["phaC"].require_pha_context is False


def test_phaA_excludes_bktb():
    families = {f.family_id: f for f in load_family_definitions()}
    assert "bktB" in families["phaA"].exclude_genes


def test_families_without_ec_skip_brenda():
    families = {f.family_id: f for f in load_family_definitions()}
    for family_id in ("phaE", "phaR_synthase", "phaG", "phaZ", "phaY", "phaP", "phaF", "phaI", "phaM", "phaD", "phaR_regulator", "phaQ"):
        assert families[family_id].ec_brenda is None, family_id


def test_monomer_and_polymer_vocab_load():
    monomers = load_monomer_vocab()
    polymers = load_polymer_vocab()
    monomer_codes = {m.code for m in monomers}
    polymer_codes = {p.code for p in polymers}
    assert {"3HB", "4HB", "3HV", "3HHx"}.issubset(monomer_codes)
    assert {"PHB", "PHBV", "P(3HB-co-4HB)"}.issubset(polymer_codes)


def test_config_file_actually_used_matches_default_path():
    assert DEFAULT_FAMILY_CONFIG_PATH.exists()
    assert DEFAULT_FAMILY_CONFIG_PATH.name == "family_definitions.yaml"
