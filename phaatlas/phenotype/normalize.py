"""Conservative monomer/polymer phenotype normalization.

The crucial evidence rule from the spec: normalization NEVER upgrades a
generic annotation into a demonstrated phenotype. Which
phenotype_evidence_level a piece of text gets is decided entirely by WHERE
it came from (source_field), not by re-reading the text itself -- e.g. a
BRENDA natural_substrates_products entry is experimental_literature by
construction (BRENDA only populates that field when the literature says so),
while a UniProt generic Rhea catalytic-activity comment is
generic_reaction_annotation even if it happens to name a specific monomer.
Only 'experimental_literature' and 'curated_specific_function' are allowed
to populate demonstrated_monomers/demonstrated_polymers; the other two
levels always leave those two fields empty, on purpose.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from phaatlas.config_loader import VocabTerm, load_monomer_vocab, load_polymer_vocab

EXPERIMENTAL_LITERATURE = "experimental_literature"
CURATED_SPECIFIC_FUNCTION = "curated_specific_function"
GENERIC_REACTION_ANNOTATION = "generic_reaction_annotation"
INFERRED_ANNOTATION = "inferred_annotation"

DEMONSTRATED_ALLOWED_LEVELS = {EXPERIMENTAL_LITERATURE, CURATED_SPECIFIC_FUNCTION}

# BRENDA fields whose entries are experimental/literature-linked by
# construction (BRENDA only populates them when a cited paper says so).
BRENDA_EXPERIMENTAL_FIELDS = {
    "natural_substrates_products",
    "substrates_products",
    "general_information",
}

# Everything else that isn't in BRENDA_EXPERIMENTAL_FIELDS is treated as
# generic/reaction-level BRENDA content (e.g. `reaction`) unless the caller
# says otherwise.


@dataclass(frozen=True)
class VocabMatcher:
    terms: tuple[VocabTerm, ...]
    _patterns: tuple[tuple[str, re.Pattern], ...]

    @staticmethod
    def build(terms: list[VocabTerm]) -> "VocabMatcher":
        patterns = []
        for term in terms:
            for syn in term.synonyms:
                # Whole-phrase, case-insensitive; escape regex metachars in
                # the synonym itself (parentheses are common, e.g. "P(3HB)").
                pat = re.compile(r"(?<![A-Za-z0-9])" + re.escape(syn) + r"(?![A-Za-z0-9])", re.IGNORECASE)
                patterns.append((term.code, pat))
        # Longest synonym first so e.g. "3-hydroxybutyrate" doesn't get
        # shadowed by a shorter unrelated match before a longer, more
        # specific one is tried.
        patterns.sort(key=lambda t: -len(t[1].pattern))
        return VocabMatcher(terms=tuple(terms), _patterns=tuple(patterns))

    def find_codes(self, text: str) -> list[str]:
        if not text:
            return []
        found: list[str] = []
        for code, pattern in self._patterns:
            if pattern.search(text) and code not in found:
                found.append(code)
        return found


_monomer_matcher: VocabMatcher | None = None
_polymer_matcher: VocabMatcher | None = None


def get_monomer_matcher() -> VocabMatcher:
    global _monomer_matcher
    if _monomer_matcher is None:
        _monomer_matcher = VocabMatcher.build(load_monomer_vocab())
    return _monomer_matcher


def get_polymer_matcher() -> VocabMatcher:
    global _polymer_matcher
    if _polymer_matcher is None:
        _polymer_matcher = VocabMatcher.build(load_polymer_vocab())
    return _polymer_matcher


@dataclass
class NormalizedPhenotype:
    phenotype_raw: str
    demonstrated_monomers: str | None
    demonstrated_polymers: str | None
    composition_text: str | None
    phenotype_evidence_level: str


def normalize_phenotype_text(raw_text: str, evidence_level: str) -> NormalizedPhenotype:
    """Extract monomer/polymer codes from raw_text using the controlled
    vocabulary, but only populate demonstrated_* if evidence_level allows it.
    The original sentence is always preserved in phenotype_raw regardless.
    """
    monomers: list[str] = []
    polymers: list[str] = []
    if evidence_level in DEMONSTRATED_ALLOWED_LEVELS and raw_text:
        monomers = get_monomer_matcher().find_codes(raw_text)
        polymers = get_polymer_matcher().find_codes(raw_text)

    composition_text = raw_text if _looks_like_composition_statement(raw_text) else None

    return NormalizedPhenotype(
        phenotype_raw=raw_text,
        demonstrated_monomers="|".join(monomers) if monomers else None,
        demonstrated_polymers="|".join(polymers) if polymers else None,
        composition_text=composition_text,
        phenotype_evidence_level=evidence_level,
    )


_COMPOSITION_HINTS = re.compile(r"\bmol ?%|\bmol\.?%|\bcomposition\b|\bmonomer\b", re.IGNORECASE)


def _looks_like_composition_statement(text: str) -> bool:
    return bool(text) and bool(_COMPOSITION_HINTS.search(text))


def brenda_field_evidence_level(field_name: str) -> str:
    """Evidence level for a raw-text statement lifted straight from a named
    BRENDA field (see BRENDA_EXPERIMENTAL_FIELDS above)."""
    if field_name in BRENDA_EXPERIMENTAL_FIELDS:
        return EXPERIMENTAL_LITERATURE
    return GENERIC_REACTION_ANNOTATION


def uniprot_function_evidence_level(function_text: str) -> str:
    """UniProt's free-text Function comment can range from protein-specific
    biology to a generic restatement of the Rhea reaction. Conservative
    heuristic: only promote to curated_specific_function when the comment
    text itself names a specific monomer/polymer from the controlled
    vocabulary (i.e. it's saying something concrete, not just "catalyzes the
    polymerization of (R)-3-hydroxyacyl-CoA monomers into PHA" in the
    abstract). Anything else stays generic_reaction_annotation -- err
    toward under-claiming per the spec's evidence rule.
    """
    if not function_text:
        return GENERIC_REACTION_ANNOTATION
    if get_monomer_matcher().find_codes(function_text) or get_polymer_matcher().find_codes(function_text):
        return CURATED_SPECIFIC_FUNCTION
    return GENERIC_REACTION_ANNOTATION
