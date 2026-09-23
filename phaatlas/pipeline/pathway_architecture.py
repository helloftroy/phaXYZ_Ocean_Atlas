"""Genome x PHA-family presence/count matrix, collapsed into short
"pathway architecture" labels (e.g. "ABC", "ABCJ") -- which combinations
of PHA genes co-occur in the same genome, how common each combination is,
and which taxa carry it.

Input: one or more <family>_unique_targets_with_metadata.tsv files (the
output of pipeline/omdb_metadata.py's enrich_unique_targets, one row per
target_id x genome, joined to that genome's GTDB taxonomy and sample
location/ecosystem). This module only reads those files -- no new network
calls or database access.

Scope note: a genome's architecture label only reflects the families
whose <family>_unique_targets_with_metadata.tsv has actually been built
(via `pha-reference omdb-enrich-metadata`) and fed into this module. A
family missing from the input set is indistinguishable here from a
genome genuinely lacking that gene -- both read as "absent". Common/rare
frequencies and "which taxa carry them" are only meaningful once most/all
families of interest have been enriched and included; with only one
family enriched so far, every genome's architecture is trivially just
that one family's own code, by construction, not a real finding.

This module does NOT infer gene order/operon structure (e.g. the
"C1-Z-C2" kind of label the pipeline may eventually want) -- only
presence/absence + counts per family per genome. Gene order would need
each hit's scaffold and within-scaffold gene index, which IS present in
the underlying NR100 cluster member IDs (<GENOME>-scaffold_<N>_<gene#>,
see pipeline/omdb_metadata.py) but is not threaded through to
unique_targets_with_metadata.tsv today -- a deliberately separate,
later piece of work.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from phaatlas.config_loader import DEFAULT_FAMILY_CONFIG_PATH, load_family_definitions

TAXONOMY_COLUMNS = [
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species",
]
LOCATION_COLUMNS = [
    "sample_id", "study_id", "latitude_degN", "longitude_degE",
    "ecosystem_type", "ecosystem_name", "ecosystem_compartment", "sample_source",
    # Only present in a *_with_metadata_depth.tsv (see pipeline/ncbi_depth.py)
    # -- a plain *_with_metadata.tsv simply has none of these three columns
    # in its header, so row.get(...) below reads them as "" gracefully.
    "depth_raw", "depth_m", "depth_zone",
]


def short_code(family_id: str) -> str:
    """phaA -> A, phaR_regulator -> RReg, phaR_synthase -> RSyn. Derived
    from family_id rather than a hardcoded table so it can never drift
    from config/family_definitions.yaml."""
    name = family_id.removeprefix("pha")
    if "_" in name:
        base, suffix = name.split("_", 1)
        return base + suffix[:3].capitalize()
    return name


def default_family_order(config_path: Path = DEFAULT_FAMILY_CONFIG_PATH) -> list[str]:
    """Canonical family ordering for architecture labels -- config file's
    own order, so labels are deterministic and match this project's
    existing family_id ordering elsewhere (e.g. `status` command output)."""
    return [f.family_id for f in load_family_definitions(config_path)]


def discover_metadata_files(results_dir: Path) -> list[Path]:
    return sorted(results_dir.glob("*_unique_targets_with_metadata.tsv"))


# Full reference-query audit, 2026-09-22: every one of the ~9,600 distinct
# reference sequences that ever recruited a target across all 15 PHA gene
# families was individually resolved against UniProt and checked against that
# family's expected function -- not a top-N sample. Two dominant, systemic
# collision sources account for most of it: (1) a real, unrelated bacterial
# multi-subunit K+/Na+-H+ antiporter locus historically also called "Pha" (for
# pH adaptation), whose subunit letters A-G collide one-for-one with this
# project's phaA-phaG gene names (confirmed explicitly in several baits' own
# UniProt names, e.g. "Pha system subunit A/B"); (2) a recurring fatty-acid-
# synthesis/phenylacetate-catabolism contaminant cluster appearing at nearly
# identical volume across phaB/D/E/F/G, suggesting one mislabeled source-genome
# cluster contaminated multiple families at once. See figures/PHAC_VALIDATION_
# FINDINGS.md section 10 (phaC) and figures/PHA_ALL_FAMILIES_REFERENCE_AUDIT.md
# (all other families) for full methodology, per-reference evidence, and the
# source audit tables. Replaces an earlier, explicitly non-exhaustive "top 5
# per family" pass.
FAMILY_BAD_QUERIES: dict[str, frozenset[str]] = {
    "phaA": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A024HB71", "A0A0A1W7Z0", "A0A0C9N453", "A0A0G3BKZ6", "A0A0H2Z7R2", "A0A0H3LQ83", "A0A0H3LV70",
        "A0A157R548", "A0A157SW25", "A0A1B0ZMR6", "A0A1L3LRE3", "A0A1R4IZE9", "A0A1S6R513", "A0A1W6NY76",
        "A0A291P4C1", "A0A2C9EF17", "A0A348FW76", "A0A380T6R6", "A0A383S180", "A0A4Y4CUL5", "A0A5B6WE42",
        "A0A640VQS7", "A0A7R6P9C9", "A0A829Y8Q7", "A0A830G059", "A0A8H9IEV8", "A0A8J2ZJ61", "A0A8J3GHQ2",
        "A0A8J3MF21", "A0A916RZH8", "A0A916TVE6", "A0A916W1I4", "A0A917ANY9", "A0A917QDP1", "A0A917SP00",
        "A0A917YI81", "A0A918IUD6", "A0A918JF96", "A0A918TGG3", "A0A918XTN1", "A0A9W6P633", "A0AA37WZ12",
        "A0AA48KIQ0", "A0AAC9UIA3", "A0AAD3NZR5", "A0AAN2PER7", "A0AAV6M9T9", "A0AB34FUA7", "A0AB34G8B8",
        "A0ABM8C0B3", "A0ABN5C9N1", "A0ABQ0Z0E9", "A0ABQ1IZU8", "A0ABQ1KEN3", "A0ABQ1QRT9", "A0ABQ3FH34",
        "A0ABQ5VKW3", "A0ABQ5ZGF2", "A0ABR9DZ94", "A0ABR9E7V5", "A0ABR9EPR7", "A0ABY1U0N4", "A0ACA8DYK6",
        "A0ACE9VIN4", "A0ACM9C487", "A8LTG8", "A9FZG0", "A9IZT8", "B0RMX9", "B9K0M4", "C4XJV8", "C6AB02",
        "C6L672", "C7BNH5", "D3BII4", "F4PPP9", "G4QHR9", "H6RQN0", "I0I5Q4", "I1DXP6", "I3TWR8", "I4EWX7",
        "I4F195", "M1PES5", "N6UDV2", "P26494", "Q162R0", "Q3IZ59", "Q54GI9", "Q5LPF5", "Q5P8E7", "Q8U5J7",
    )),
    "phaB": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A024HBH7", "A0A0A1W7Z0", "A0A0C4WQB3", "A0A0C9N453", "A0A0H2UXA9", "A0A0Q0EQ65", "A0A0W0RZ30",
        "A0A0W0S9K9", "A0A0W0SBJ7", "A0A0W0SUA3", "A0A0W0SXN2", "A0A0W0VCH5", "A0A0W0VIF6", "A0A0W0W157",
        "A0A0W0WP20", "A0A0W0Y892", "A0A0W0Z7I0", "A0A0W0ZKW2", "A0A0W0ZSY3", "A0A1L3LJX0", "A0A1L7LHN0",
        "A0A1M7E650", "A0A1S6R4X7", "A0A239SYK5", "A0A2C9EF77", "A0A378I622", "A0A378IEF4", "A0A378IZY2",
        "A0A378L4F3", "A0A378LM20", "A0A380L0L9", "A0A380T817", "A0A3G8JNJ5", "A0A3M2RBN0", "A0A447Z355",
        "A0A498QZ26", "A0A4Y3PIX3", "A0A4Y3PQZ5", "A0A511UU87", "A0A511UXV7", "A0A511V5H5", "A0A511VE22",
        "A0A511Z4C0", "A0A511ZDP0", "A0A5M3XDH4", "A0A5S9R3N5", "A0A6F8T9C2", "A0A6N4SV21", "A0A7D5Z473",
        "A0A7U7ENM0", "A0A7U7EP66", "A0A7U9XVC3", "A0A7Z7NB23", "A0A7Z7NCJ2", "A0A8J2ZXE8", "A0A917AUY8",
        "A0A917EEV1", "A0A917G8D4", "A0A917JVZ2", "A0A917JYU2", "A0A918DTI0", "A0A919WE21", "A0A9P3UZN1",
        "A0A9Q0N7F9", "A0AA50DMA9", "A0AAN2TRI8", "A0ABD7NCB0", "A0ABN5XY08", "A0ABQ0RX64", "A0ABQ1K9B8",
        "A0ABQ2D8X9", "A0ABQ2DIU2", "A0ABQ2NY13", "A0ABQ2ZNX6", "A0ABQ3MZK0", "A0ABQ4K975", "A0ABQ4L407",
        "A0ABQ5NJ19", "A0ABQ5PET1", "A0ABQ5TKV9", "A0ABQ5ZRP4", "A0ABS0KSJ6", "A0ABS6MTW5", "A0ABS9FLK7",
        "A0ABT4XJN4", "A0ABT6IF28", "A0ABU7N1Q0", "A0ABU8PT37", "A0ABV3YTH4", "A0ABV6SQF6", "A0ABW8WA51",
        "A0ABX8YM61", "A0ABY5EMB4", "A3CQ62", "A4FE21", "A4VHE0", "C9JRZ8", "D3HNN7", "D5DZ99", "E0PBP1",
        "H6RQN1", "I4F194", "K0CJG3", "P71534", "P9WGT3", "Q0C0M1", "Q0VS74", "Q52978", "Q54Q31", "Q5M5S3",
        "Q87UZ8", "Q8DR19", "Q8KR81", "Q9Z3Y5", "T1ZEW6",
    )),
    "phaC": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0D6JHM2", "A0A0G3BRR5", "A0A0H2Z7T5", "A0A0H3LUG2", "A0A0H3LVF9", "A0A0R0DVU2", "A0A0R4J7U1",
        "A0A157R4S9", "A0A157SVZ8", "A0A1B0ZMT2", "A0A1E8EW93", "A0A1G5SDS9", "A0A1S6R531", "A0A1W6NY09",
        "A0A348FW77", "A0A480AP87", "A0A6S6XVD8", "A0A7U7GEB2", "A0A829YA03", "A0A8J3E0L9", "A0A8S2BDR2",
        "A0A916RII9", "A0A916RIM5", "A0A916S023", "A0A916ZXM0", "A0A917QD60", "A0A918ITM8", "A0A918W8J3",
        "A0ABN5C053", "A0ABN5CAC2", "A0ABN6T5N2", "A0ABP2DAC8", "A0ABQ1J139", "A0ABQ1KAJ0", "A0ABQ2CQ11",
        "A0ABQ2WTK0", "A0ABQ5TXW2", "A0ABQ5VKV7", "A9CK78", "A9FZG2", "A9IP81", "B2FKM3", "B6JCB5",
        "C7BNH2", "C8YNX0", "I3R9Z3", "I3TJ70", "M0I0U4", "M0ICI7", "N6UCP8", "N6VCC7", "P45367", "P45372",
        "P73389", "Q162R1", "Q17V39", "Q17V43", "Q2KZX2", "Q4KDG6", "Q52980", "Q5P8E8", "Q5UYM1", "Q6FF34",
        "Q6RI97", "Q8KR79", "Q92RA1", "Q9F5P9",
        # Added 2026-09-22: full direct InterPro re-verification of every previously
        # ON_TARGET-classified phaC reference (1,867 queries), after the original
        # 8-query "FUSION-ABC_transporter" exclusion turned out to be non-exhaustive.
        # Root cause: interpro_full_reaudit.py's GOOD_DOMAINS never had a 'phaC' entry
        # (every other family did), so phaC references were never checked for the
        # actual confirmed marker domain (IPR051321) -- only screened against a
        # negative blacklist, which a differently-wrong protein passes easily. Checked
        # directly against live UniProt/InterPro data, not just recruitment statistics;
        # each failure mode spot-verified (misannotated non-PhaC enzymes with no
        # relation to IPR051321, e.g. acyl-CoA synthetases and polysaccharide-
        # biosynthesis proteins; multi-domain fusions carrying a real but minority
        # PhaC domain alongside an unrelated dominant one; implausibly short/fragment
        # or organism-mismatched entries with no InterPro annotation at all). See
        # PHA_CLEAN_RESULTS.md for the full verification writeup.
        "A0A090R7C3", "A0A097ERZ2", "A0A0C3DJ16", "A0A0D8PXE4", "A0A0G3WHX8", "A0A0J1GRI2", "A0A0J1K4W4", "A0A0T7DNW4", "A0A150WMA4",
        "A0A150WNK1", "A0A151AEV0", "A0A151AIL1", "A0A151L0P1", "A0A178JFB0", "A0A1B9QV15", "A0A1G8Y6C9", "A0A1L3GI38", "A0A1Q9G762",
        "A0A1T5FPD3", "A0A1X7AFX8", "A0A1Y4DIX0", "A0A1Y6KRP4", "A0A1Y6M9P7", "A0A240EG07", "A0A2M8H6D7", "A0A2N4UXW8", "A0A2P8GDY2",
        "A0A2S7V7Z2", "A0A2T3NAI5", "A0A2T3PRJ9", "A0A330LTS8", "A0A345S6Y8", "A0A3G8JFS6", "A0A3Q9K0Y4", "A0A3S9XZH6", "A0A444JWC5",
        "A0A498Q0R0", "A0A4D6LUZ9", "A0A4R1K402", "A0A4R2NZ35", "A0A4U3FLS4", "A0A4Y3HY53", "A0A4Y3IMB5", "A0A4Y3IRN2", "A0A4Y8WJ87",
        "A0A556QKC2", "A0A557ST61", "A0A5J6WQD8", "A0A5J6WST3", "A0A5J6WUA4", "A0A5Q0TAN5", "A0A5S9QVF8", "A0A6A7GGQ7", "A0A6P1M6T9",
        "A0A7J5AGP2", "A0A7X1E4U9", "A0A7Z2T3G5", "A0A7Z7ILZ6", "A0A812VQC3", "A0A813BGK3", "A0A853R7Z3", "A0A8S2BLM7", "A0A918JWF9",
        "A0A918JXY7", "A0A9W6GI39", "A0A9W6GL72", "A0A9W6GNE3", "A0A9W6LLT3", "A0A9X1WD73", "A0A9X3HW65", "A0AAD1BY69", "A0AAI7ZCY1",
        "A0AAN1CWM3", "A0AAN4VX63", "A0AAU8BIW8", "A0AAU9D2C6", "A0AAV5NSL7", "A0AAX0Z1X0", "A0AAX2LJW2", "A0ABM7GZ97", "A0ABM8ZQR0",
        "A0ABN5JGB5", "A0ABN8E4V9", "A0ABN8EHT7", "A0ABN8JP91", "A0ABP7VE29", "A0ABP8QGS0", "A0ABP9S7A1", "A0ABQ1I0H1", "A0ABQ5Y5A3",
        "A0ABQ6DXS1", "A0ABR7J2E5", "A0ABS2G0X7", "A0ABS2HP62", "A0ABT1N2I8", "A0ABT7QAT8", "A0ABU4W9W4", "A0ABU7G4P0", "A0ABU9FSZ5",
        "A0ABU9HU42", "A0ABU9JAJ3", "A0ABW9GBA2", "A0ABX0DA28", "A0ABX1U6D0", "A0ABX3AZ03", "A0ABX4XAX1", "A0ABX5GGB4", "A0ABX5H8K4",
        "A0ABX5H8M3", "A0ABX9KEK0", "A0ABX9KJH0", "A0ABX9KJK0", "A0ABY1HDR6", "A0ABY3S282", "A0ABY5G6F1", "A0ABZ0F622", "A0ACE0BGR7",
        "A0ACH1R8B4", "A0ACH1RK58", "A0ACH1RMQ7", "A0ACH3UPZ5", "A0ACM8RTF1", "A0ACM8RWN5", "A0ACN3LVC5", "A0ACN3P6Y2", "A0ACN3VS62",
        "A0ACN4YJT2", "D0YZB4", "E1SWP3", "E3HAX3", "E8M187", "K1KB67", "K6DJS5", "M0DCW3", "N9VHT5",
        "P76108", "Q3A8F5", "Q5E3C0", "Q87KX4", "Q8GFF3", "Q8GI83", "Q8PD95", "Q9KNR9", "R1F852",
        "U3B5V2", "U4KH66", "U5NXH1", "U7D663", "U7VCC2", "V5FG99", "V5FMS0", "W5U290",
    )),
    "phaD": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0A1W8R3", "A0A0B3RMS8", "A0A0C9NDJ2", "A0A0D5YLM2", "A0A0G3BU60", "A0A0H2Z7H6", "A0A0H2ZHH5",
        "A0A0H2ZIQ8", "A0A0H3LNR5", "A0A0H3LYS4", "A0A0H3M462", "A0A0S1B3G1", "A0A157R312", "A0A157SWC2",
        "A0A1B0ZMX7", "A0A1G5SDU2", "A0A1L3LR96", "A0A1S6R518", "A0A1U9JSD8", "A0A238D4P8", "A0A291P4C2",
        "A0A2I0CNH5", "A0A2I0CNJ6", "A0A380T200", "A0A383RX56", "A0A3M4VJL7", "A0A4Y3WK08", "A0A511D6A9",
        "A0A511DIR1", "A0A512AL92", "A0A512DZB6", "A0A640VU20", "A0A6J4ED71", "A0A6S6XS24", "A0A7I7L295",
        "A0A7I7W9X9", "A0A7I7YWX1", "A0A7I9VBL0", "A0A7R6PK29", "A0A7U7GAQ3", "A0A7Z7HP56", "A0A822V230",
        "A0A829Y8R0", "A0A846LRY4", "A0A8H9IG67", "A0A8I0MZB8", "A0A8J2YMT3", "A0A8J3AYK5", "A0A8J3CQW1",
        "A0A8J3DLT9", "A0A8J3EFA8", "A0A8J3HCA4", "A0A8S2BK31", "A0A916BB52", "A0A916RJF9", "A0A916RYU2",
        "A0A916TP91", "A0A916Y8Q2", "A0A916ZX56", "A0A917BP09", "A0A917EMW1", "A0A917QEP1", "A0A917SMP6",
        "A0A917YJR0", "A0A918DMN1", "A0A918IT73", "A0A918JGR6", "A0A918N0U3", "A0A918VKV3", "A0A918W9Q2",
        "A0A918WHT0", "A0A918XTV9", "A0A919KIC8", "A0A921NT84", "A0A9P1JVC8", "A0A9W6J7R7", "A0A9W6K1I0",
        "A0A9W6K9Z6", "A0A9W6KCF7", "A0A9W6L3I1", "A0AA37X0H6", "A0AA48KMR1", "A0AA48KU51", "A0AA86GJK5",
        "A0AAC9UH89", "A0AAD3NZL4", "A0AAI8EP32", "A0AAQ1P6Q0", "A0AAQ2ITT7", "A0ABM7IGG8", "A0ABM7U6T1",
        "A0ABM8HBV5", "A0ABM9DVQ9", "A0ABQ0Z0L1", "A0ABQ1F8H7", "A0ABQ1HPT4", "A0ABQ1J1D3", "A0ABQ1JR75",
        "A0ABQ1KAE6", "A0ABQ1LHW5", "A0ABQ1PFU8", "A0ABQ1QU33", "A0ABQ1RME2", "A0ABQ2EB88", "A0ABQ2KEQ3",
        "A0ABQ2WPX8", "A0ABQ3CZ04", "A0ABQ3FH65", "A0ABQ3HQ89", "A0ABQ3J2P5", "A0ABQ4PXY8", "A0ABQ4PY40",
        "A0ABQ5PEP1", "A0ABQ5UKG7", "A0ABQ5VBH9", "A0ABQ5VL13", "A0ABQ5ZFM1", "A0ABQ6AN02", "A0ABR9DZ92",
        "A0ABR9E7U9", "A0ABR9EPS0", "A0ABY1U5B8", "A0ACA8DYD4", "A0ACD7VSA9", "A0ACE3TUV5", "A0ACE9VIY1",
        "A0ACN7GMZ3", "A8LTH0", "A9FZG5", "A9IP78", "A9IZT4", "A9WT02", "B0RMX7", "B6JCB4", "B8KSS7",
        "B9K0M2", "D5AV78", "D5P5L0", "D6CSQ1", "F0KLB9", "G4QHR7", "G7ZDE2", "H5TAX5", "I0HVN1", "I1DXP8",
        "I2GCQ8", "K0C4S0", "M1N5G6", "N6UET0", "N6UJG7", "N6VI43", "Q0C0L9", "Q0VTK7", "Q162R2", "Q2KZX3",
        "Q3IZ57", "Q4KDG7", "Q4KJK8", "Q52981", "Q5LPF3", "Q5P8E9", "Q6FF35", "Q6LR00", "Q82TI4", "Q87UZ6",
        "S0FU87", "U2ENQ3", "U2YMZ3", "U2ZQ88", "V6CKT9", "W6M4I4",
    )),
    "phaE": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0A1W8C2", "A0A0C9MUK7", "A0A0D5YN30", "A0A0G3BKT7", "A0A0H2Z803", "A0A0H3LNR3", "A0A0K0GH14",
        "A0A157R568", "A0A157SVZ7", "A0A1B0ZMR0", "A0A1S6R4Y9", "A0A1U9JSE9", "A0A1W6NXY9", "A0A291P496",
        "A0A2K4X7N1", "A0A380T3U9", "A0A383RVH1", "A0A480APD1", "A0A512C0Z7", "A0A512DTH7", "A0A512JBG9",
        "A0A5M3PVT1", "A0A640VT22", "A0A6S6XXF2", "A0A7Z7MUP1", "A0A822V3U9", "A0A829Y950", "A0A8H9IKX3",
        "A0A8I0T561", "A0A8J2YMM8", "A0A8J3DLJ4", "A0A8J3EGN6", "A0A8J3GLK8", "A0A8J3HCA5", "A0A916RIC3",
        "A0A916S1X6", "A0A916TUK9", "A0A916Y9R0", "A0A917AN34", "A0A917F6I5", "A0A917SPS8", "A0A917V771",
        "A0A917YJ66", "A0A918ITR4", "A0A918THA9", "A0A918WAM1", "A0A918XV40", "A0A919F8P5", "A0A921NV37",
        "A0A9P1JV73", "A0A9W6K077", "A0AA37X264", "A0AA86GK14", "A0AAC9UE42", "A0AAD3NZU3", "A0AAI8ERD6",
        "A0AAN2PEH8", "A0ABN5C6I4", "A0ABN5C9T5", "A0ABN6P6G3", "A0ABN8JS05", "A0ABQ1HPM0", "A0ABQ1IY91",
        "A0ABQ1KCN1", "A0ABQ1QUM1", "A0ABQ2CPW5", "A0ABQ2E6Q1", "A0ABQ2YLM3", "A0ABQ2YZE9", "A0ABQ3FAQ6",
        "A0ABQ3FHI1", "A0ABQ3KZC5", "A0ABQ5UJ05", "A0ABQ5VKT9", "A0ABQ5ZJL6", "A0ABR9DZ93", "A0ABR9E7U8",
        "A0ACA8DXD1", "A0ACD7VRJ6", "A0ACE2TPR6", "A1K857", "A8LTH1", "A9FZG8", "A9WT01", "B0RMX6",
        "B6JCB3", "B9K0M1", "C3M903", "D5AV77", "D5RP21", "F0KLB8", "G4QHR6", "G7ZDE3", "H5TAX6", "H8FU32",
        "I0HVN2", "I1DXP9", "K0C9A7", "K6Y558", "M7NK86", "O33469", "Q162R3", "Q2KZX6", "Q3IZ56", "Q4KDG8",
        "Q52982", "Q5LPF2", "Q5P8F0", "Q6FF36", "Q7D1B4", "Q9KGB3", "S0G3B7", "U2Z794", "U2ZXN4",
    )),
    "phaF": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0A1W9I8", "A0A0C9M3G6", "A0A0G3BQ06", "A0A0H2Z6P3", "A0A0H3G746", "A0A0H3LWU6", "A0A0K0GGZ3",
        "A0A0R0DWS0", "A0A157R5B2", "A0A1A8XVB2", "A0A1B0ZMQ3", "A0A1G5SEE6", "A0A1L3LRB2", "A0A1S6R4W7",
        "A0A1S7PW04", "A0A1S7TMU0", "A0A1U9JSF0", "A0A1W6NYL6", "A0A291P4A3", "A0A2C9ELD4", "A0A380T443",
        "A0A512C1F4", "A0A512DTJ3", "A0A512DZZ4", "A0A512JBI3", "A0A640VRC7", "A0A6S6Y8H9", "A0A6V7MKE0",
        "A0A7R6P258", "A0A7U7GAE7", "A0A822UYT9", "A0A8J2XV09", "A0A8J3EFY4", "A0A8J3GLQ7", "A0A8J3HBZ0",
        "A0A8S2BM45", "A0A916BB50", "A0A916TU67", "A0A916W1M2", "A0A916W968", "A0A916Y8B0", "A0A916ZJ55",
        "A0A917AQA8", "A0A917E9U7", "A0A917F435", "A0A917SM89", "A0A917V765", "A0A917YKC2", "A0A918MJQ9",
        "A0A918T2S0", "A0A918TFZ7", "A0A919CQG7", "A0A919F8H1", "A0A921NRE0", "A0A9P1JV81", "A0A9W5B2L2",
        "A0A9W6K3W3", "A0AA37TQG3", "A0AA48H9J5", "A0AA86GJM9", "A0AAC9UIG9", "A0AAD3NZZ4", "A0ABM9SEM5",
        "A0ABN4AYB3", "A0ABN5C3F8", "A0ABN8JTP2", "A0ABP2BCP7", "A0ABP2BSI5", "A0ABQ0Z0K1", "A0ABQ1HNJ8",
        "A0ABQ1J159", "A0ABQ1K9K1", "A0ABQ1QTH3", "A0ABQ1RGW9", "A0ABQ2E6T1", "A0ABQ3FH10", "A0ABQ5UJ33",
        "A0ABQ5VKT0", "A0ABQ5ZDU9", "A0ABQ5ZWS1", "A0ABR9E7V6", "A0ACD7VS64", "A0ACE3TV65", "A0ACE9VIQ9",
        "A4BLD0", "A8LTH2", "A9CK79", "A9FZH2", "A9WT00", "B0RMX5", "B6JCB2", "B9K0M0", "C3M904", "D5RP22",
        "F0KLB7", "F9Y398", "G4QHR5", "G7ZDE4", "H5TAX7", "I0HVN3", "I1DXQ0", "K0PY53", "K6ZG59", "N6VAV8",
        "O68038", "Q162R4", "Q2KZX8", "Q3IZ55", "Q52983", "Q5LPF1", "Q5P8F1", "Q6FF37", "Q7AK13", "S0G287",
        "U2ZRA8", "U3A5Y7", "V5VHB1", "W6M9Y3",
    )),
    "phaG": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0A1W7Z5", "A0A0C9N447", "A0A0H2Z7R7", "A0A0R0DL19", "A0A0R4J8F1", "A0A0T7CPX3", "A0A157R4V6",
        "A0A157SWC8", "A0A1A8XXC0", "A0A1B0ZMZ1", "A0A1L3LR79", "A0A1S6R4Y2", "A0A1U9JSD6", "A0A1W6NXZ7",
        "A0A291P4D0", "A0A2K4X7P8", "A0A380T1W9", "A0A383RWI4", "A0A512DTJ4", "A0A512JBP1", "A0A7R6SSP4",
        "A0A829YAN7", "A0A8J2ULC6", "A0A8J2YNA7", "A0A8J2ZIY6", "A0A8J3AT24", "A0A8J3DKA5", "A0A8J3GMR3",
        "A0A8J3H8U5", "A0A916F8D0", "A0A916WKI8", "A0A916Y9H6", "A0A917WAU2", "A0A917YJ11", "A0A918JF92",
        "A0A918MJB1", "A0A918TGH0", "A0A921TDB4", "A0A9P1JVG7", "A0A9W5B1E8", "A0A9W6MYD3", "A0AA37TU55",
        "A0AA48H2L1", "A0AA86GKT9", "A0AAC9UH03", "A0AAD3NZP3", "A0AAN2PEW9", "A0AAQ1P865", "A0AAQ2EY36",
        "A0ABM8C0C1", "A0ABM8FCC1", "A0ABQ1KDC8", "A0ABQ1KUU4", "A0ABQ1QV65", "A0ABQ2WS28", "A0ABQ3FH03",
        "A0ABQ3KZ18", "A0ABQ5ZJA7", "A0ABR9DZA8", "A0ABR9E7V8", "A0ACA8DX70", "A0ACD7VS75", "A0ACE3TVA2",
        "A0ACE9VID8", "A8LTH3", "A9FZH5", "A9WSZ9", "B0RMX4", "B6JCB1", "B9K0L9", "C3M905", "D5AV75",
        "F0KLB6", "G4QHR4", "G7ZCL7", "I0HVN4", "I1DXQ1", "M1PER9", "Q2KZX9", "Q3IZ54", "Q4KDH0", "Q5LPF0",
        "Q5P8F2", "Q6FF38", "Q9KGB1", "Q9Z3Q3", "S0FT44", "U2YNK4", "U3B9V8",
    )),
    "phaI": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A1S6R537",
    )),
    "phaJ": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A1S6R502", "A0A2S2DUT1", "A0A369QCL2", "F0EWD4", "F2B9Z4", "G4CK88", "G4CPY8", "Q7NZL6",
        "Q9LBK2",
    )),
    "phaR_regulator": frozenset(f"UNIPROT:{acc}" for acc in (
        "Q1JR79", "Q8KRE8",
    )),
    "phaR_synthase": frozenset(f"UNIPROT:{acc}" for acc in (
        "Q8KRE8",
    )),
    "phaZ": frozenset(f"UNIPROT:{acc}" for acc in (
        "A0A0F3KLG9",
    )),
}

# backwards-compat alias -- cli.py's cluster-ecology/cluster-map-points commands
# import this name directly for phaC specifically
PHAC_BAD_QUERIES = FAMILY_BAD_QUERIES["phaC"]


@dataclass
class GenomeRecord:
    genome: str
    taxonomy: dict[str, str] = field(default_factory=dict)
    location: dict[str, str] = field(default_factory=dict)
    family_target_ids: dict[str, set[str]] = field(default_factory=dict)  # family_id -> distinct target_ids hitting this genome

    def family_counts(self) -> dict[str, int]:
        return {fam: len(ids) for fam, ids in self.family_target_ids.items()}


def load_genome_records(metadata_paths: list[Path]) -> dict[str, GenomeRecord]:
    """Reads every given <family>_unique_targets_with_metadata.tsv and
    aggregates per-genome family hit counts + taxonomy/location (taken
    from whichever row is seen first for that genome -- these are
    genome/sample facts, not family-specific, so any row agrees)."""
    genomes: dict[str, GenomeRecord] = {}
    for path in metadata_paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome = row.get("genome", "")
                if not genome:
                    continue  # no genome resolved for this target_id -- nothing to attribute
                family_id = row["pha_family"]
                target_id = row["target_id"]
                if row.get("best_query") in FAMILY_BAD_QUERIES.get(family_id, ()):
                    continue

                rec = genomes.get(genome)
                if rec is None:
                    rec = GenomeRecord(
                        genome=genome,
                        taxonomy={c: row.get(c, "") for c in TAXONOMY_COLUMNS},
                        location={c: row.get(c, "") for c in LOCATION_COLUMNS},
                    )
                    genomes[genome] = rec
                rec.family_target_ids.setdefault(family_id, set()).add(target_id)
    return genomes


def architecture_label(counts: dict[str, int], family_order: list[str]) -> str:
    return "".join(short_code(fam) for fam in family_order if counts.get(fam, 0) > 0)


def write_genome_family_matrix(
    genomes: dict[str, GenomeRecord],
    family_order: list[str],
    out_path: Path,
) -> int:
    """One row per genome: taxonomy/location + one count column per
    family in family_order + the collapsed architecture label. Returns
    the number of rows written."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        ["genome"] + TAXONOMY_COLUMNS + LOCATION_COLUMNS
        + [f"n_{fam}" for fam in family_order] + ["n_families_present", "architecture"]
    )
    n = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for genome in sorted(genomes):
            rec = genomes[genome]
            counts = rec.family_counts()
            row = {"genome": genome, **rec.taxonomy, **rec.location}
            for fam in family_order:
                row[f"n_{fam}"] = counts.get(fam, 0)
            row["n_families_present"] = sum(1 for fam in family_order if counts.get(fam, 0) > 0)
            row["architecture"] = architecture_label(counts, family_order)
            writer.writerow(row)
            n += 1
    return n


@dataclass
class ArchitectureStats:
    architecture: str
    families_present: list[str]
    n_genomes: int
    pct_of_genomes: float
    top_genera: list[tuple[str, int]]
    top_species: list[tuple[str, int]]
    n_distinct_studies: int
    top_studies: list[tuple[str, int]]


def summarize_architectures(
    genomes: dict[str, GenomeRecord],
    family_order: list[str],
    top_n: int = 5,
) -> list[ArchitectureStats]:
    """One ArchitectureStats per distinct architecture found, sorted most
    to least common (ties broken alphabetically for determinism)."""
    total = len(genomes)
    by_arch: dict[str, list[GenomeRecord]] = {}
    for rec in genomes.values():
        arch = architecture_label(rec.family_counts(), family_order)
        by_arch.setdefault(arch, []).append(rec)

    results: list[ArchitectureStats] = []
    for arch, recs in by_arch.items():
        genera = Counter(r.taxonomy.get("gtdb_genus", "") for r in recs if r.taxonomy.get("gtdb_genus"))
        species = Counter(r.taxonomy.get("gtdb_species", "") for r in recs if r.taxonomy.get("gtdb_species"))
        studies = Counter(r.location.get("study_id", "") for r in recs if r.location.get("study_id"))
        families_present = [fam for fam in family_order if short_code(fam) in _split_codes(arch, family_order)]
        results.append(ArchitectureStats(
            architecture=arch,
            families_present=families_present,
            n_genomes=len(recs),
            pct_of_genomes=100.0 * len(recs) / total if total else 0.0,
            top_genera=genera.most_common(top_n),
            top_species=species.most_common(top_n),
            n_distinct_studies=len(studies),
            top_studies=studies.most_common(top_n),
        ))
    results.sort(key=lambda s: (-s.n_genomes, s.architecture))
    return results


def _split_codes(arch: str, family_order: list[str]) -> set[str]:
    """Architecture labels are a straight concatenation of variable-length
    short codes (e.g. RReg is 4 chars) -- reconstructing which families
    are present from the label alone means matching against the known
    code set greedily longest-first, rather than assuming a fixed width."""
    codes_by_length = sorted({short_code(f) for f in family_order}, key=len, reverse=True)
    present = set()
    i = 0
    while i < len(arch):
        for code in codes_by_length:
            if arch.startswith(code, i):
                present.add(code)
                i += len(code)
                break
        else:
            raise ValueError(f"Could not decompose architecture label {arch!r} using known family codes")
    return present


def write_architecture_summary(stats: list[ArchitectureStats], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "architecture", "families_present", "n_genomes", "pct_of_genomes",
            "n_distinct_studies", "top_genera", "top_species", "top_studies",
        ])
        for s in stats:
            writer.writerow([
                s.architecture,
                ",".join(s.families_present),
                s.n_genomes,
                f"{s.pct_of_genomes:.2f}",
                s.n_distinct_studies,
                "; ".join(f"{g} ({n})" for g, n in s.top_genera),
                "; ".join(f"{sp} ({n})" for sp, n in s.top_species),
                "; ".join(f"{st} ({n})" for st, n in s.top_studies),
            ])
