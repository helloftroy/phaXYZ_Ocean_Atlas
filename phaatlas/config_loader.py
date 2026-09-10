"""Loads config/family_definitions.yaml and config/vocab/*.yaml.

This is the single place that turns the YAML family config into Python
objects -- phaatlas's sources/pipeline modules never read family
specifics (EC numbers, gene aliases, search terms) from anywhere else, so
adding or editing a family is a config-only change.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FAMILY_CONFIG_PATH = REPO_ROOT / "config" / "family_definitions.yaml"
DEFAULT_MONOMER_VOCAB_PATH = REPO_ROOT / "config" / "vocab" / "monomers.yaml"
DEFAULT_POLYMER_VOCAB_PATH = REPO_ROOT / "config" / "vocab" / "polymers.yaml"


@dataclass(frozen=True)
class FamilyDefinition:
    family_id: str
    display_name: str
    role: str
    genes: tuple[str, ...]
    exclude_genes: tuple[str, ...]
    protein_name_terms: tuple[str, ...]
    ec_brenda: str | None
    require_pha_context: bool
    disambiguation_terms: tuple[str, ...]
    broad_uniprot_search: bool
    raw: dict = field(compare=False)

    def raw_json(self) -> str:
        return json.dumps(self.raw, sort_keys=True)


def load_family_definitions(path: Path = DEFAULT_FAMILY_CONFIG_PATH) -> list[FamilyDefinition]:
    with open(path, "r") as f:
        doc = yaml.safe_load(f)

    families = []
    seen_ids: set[str] = set()
    for entry in doc.get("families", []):
        family_id = entry["family_id"]
        if family_id in seen_ids:
            raise ValueError(f"duplicate family_id '{family_id}' in {path}")
        seen_ids.add(family_id)
        families.append(
            FamilyDefinition(
                family_id=family_id,
                display_name=entry["display_name"],
                role=entry["role"],
                genes=tuple(entry.get("genes") or []),
                exclude_genes=tuple(entry.get("exclude_genes") or []),
                protein_name_terms=tuple(entry.get("protein_name_terms") or []),
                ec_brenda=entry.get("ec_brenda"),
                require_pha_context=bool(entry.get("require_pha_context", False)),
                disambiguation_terms=tuple(entry.get("disambiguation_terms") or []),
                broad_uniprot_search=bool(entry.get("broad_uniprot_search", False)),
                raw=entry,
            )
        )

    # phaR_synthase and phaR_regulator sharing gene "phaR" is expected and
    # intentional (see family_definitions.yaml); any OTHER pair of families
    # silently sharing an exact gene name is more likely a config typo, so
    # flag it instead of letting it silently double-assign proteins later.
    gene_to_families: dict[str, list[str]] = {}
    for fam in families:
        for gene in fam.genes:
            gene_to_families.setdefault(gene.lower(), []).append(fam.family_id)
    expected_shared = {frozenset({"phaR_synthase", "phaR_regulator"})}
    for gene, fam_ids in gene_to_families.items():
        if len(fam_ids) > 1 and frozenset(fam_ids) not in expected_shared:
            raise ValueError(
                f"gene '{gene}' is claimed by multiple families {fam_ids} "
                "without an expected disambiguation pair -- check family_definitions.yaml"
            )

    return families


def get_family(family_id: str, path: Path = DEFAULT_FAMILY_CONFIG_PATH) -> FamilyDefinition:
    for fam in load_family_definitions(path):
        if fam.family_id == family_id:
            return fam
    raise KeyError(f"unknown family_id '{family_id}' (not found in {path})")


@dataclass(frozen=True)
class VocabTerm:
    code: str
    synonyms: tuple[str, ...]


def load_vocab(path: Path) -> list[VocabTerm]:
    with open(path, "r") as f:
        doc = yaml.safe_load(f)
    key = "monomers" if "monomers" in doc else "polymers"
    terms = []
    for entry in doc[key]:
        terms.append(VocabTerm(code=str(entry["code"]), synonyms=tuple(entry.get("synonyms") or [])))
    return terms


def load_monomer_vocab(path: Path = DEFAULT_MONOMER_VOCAB_PATH) -> list[VocabTerm]:
    return load_vocab(path)


def load_polymer_vocab(path: Path = DEFAULT_POLYMER_VOCAB_PATH) -> list[VocabTerm]:
    return load_vocab(path)
