from phaatlas.phenotype.normalize import (
    CURATED_SPECIFIC_FUNCTION,
    EXPERIMENTAL_LITERATURE,
    GENERIC_REACTION_ANNOTATION,
    INFERRED_ANNOTATION,
    brenda_field_evidence_level,
    normalize_phenotype_text,
    uniprot_function_evidence_level,
)


def test_brenda_experimental_fields_are_classified_experimental_literature():
    assert brenda_field_evidence_level("natural_substrates_products") == EXPERIMENTAL_LITERATURE
    assert brenda_field_evidence_level("substrates_products") == EXPERIMENTAL_LITERATURE
    assert brenda_field_evidence_level("general_information") == EXPERIMENTAL_LITERATURE


def test_unlisted_brenda_field_defaults_to_generic():
    assert brenda_field_evidence_level("reaction") == GENERIC_REACTION_ANNOTATION


def test_demonstrated_fields_populated_for_experimental_literature():
    text = "accumulates poly(3-hydroxybutyrate-co-4-hydroxybutyrate) with 4-hydroxybutyrate monomer up to 87 mol%"
    result = normalize_phenotype_text(text, EXPERIMENTAL_LITERATURE)
    assert result.demonstrated_monomers is not None
    assert "4HB" in result.demonstrated_monomers.split("|")
    assert result.demonstrated_polymers == "P(3HB-co-4HB)"
    assert result.phenotype_raw == text  # original sentence always preserved
    assert result.composition_text is not None  # "mol%" hint


def test_demonstrated_fields_populated_for_curated_specific_function():
    text = "Polymerizes (R)-3-hydroxybutyryl-CoA to create polyhydroxybutyrate (PHB)"
    result = normalize_phenotype_text(text, CURATED_SPECIFIC_FUNCTION)
    assert result.demonstrated_polymers == "PHB"


def test_crucial_evidence_rule_generic_annotation_never_populates_demonstrated_fields():
    # Same text, but sourced from a generic Rhea/catalytic-activity style
    # annotation -- must NOT populate demonstrated_monomers/polymers even
    # though the text itself names a specific polymer.
    text = "polymerizes (R)-3-hydroxybutyrate into poly(3-hydroxybutyrate) (PHB)"
    generic = normalize_phenotype_text(text, GENERIC_REACTION_ANNOTATION)
    assert generic.demonstrated_monomers is None
    assert generic.demonstrated_polymers is None
    assert generic.phenotype_raw == text  # raw sentence still preserved


def test_crucial_evidence_rule_inferred_annotation_never_populates_demonstrated_fields():
    text = "poly(3-hydroxybutyrate) PHB 3-hydroxybutyrate"
    inferred = normalize_phenotype_text(text, INFERRED_ANNOTATION)
    assert inferred.demonstrated_monomers is None
    assert inferred.demonstrated_polymers is None


def test_uniprot_function_text_with_specific_monomer_is_curated_specific_function():
    text = "Polymerizes (R)-3-hydroxybutyryl-CoA to create polyhydroxybutyrate (PHB)"
    assert uniprot_function_evidence_level(text) == CURATED_SPECIFIC_FUNCTION


def test_uniprot_function_text_without_specific_monomer_is_generic():
    text = "Catalyzes the polymerization of (R)-3-hydroxyacyl-CoA monomers into PHA"
    # deliberately vague/generic -- no controlled-vocabulary monomer/polymer term
    assert uniprot_function_evidence_level(text) == GENERIC_REACTION_ANNOTATION


def test_empty_text_never_crashes_and_is_generic():
    assert uniprot_function_evidence_level("") == GENERIC_REACTION_ANNOTATION
    assert uniprot_function_evidence_level(None) == GENERIC_REACTION_ANNOTATION


def test_monomer_matching_is_whole_phrase_not_substring():
    # "PHB" should not spuriously match inside an unrelated longer token
    text = "the alphabet includes PHBX99 as an internal code, unrelated to PHB itself"
    result = normalize_phenotype_text(text, EXPERIMENTAL_LITERATURE)
    assert result.demonstrated_polymers == "PHB"  # only the real, boundary-matched occurrence counts
