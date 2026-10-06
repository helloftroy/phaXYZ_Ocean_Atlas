"""One canonical habitat label set, shared by every figure that breaks phaC
down by ocean habitat.

The bar chart (plot_phac_pct_by_ocean_habitat.py) used to carry its own
DISPLAY dict renaming OMDB's `ecosystem_compartment` values, while the
phylum-controlled test (phac_habitat_phylum_controlled_test.py) wrote the raw
values straight out. Put side by side, the same habitat appeared under two
names in two figures -- "Marine algae thallus" against "Algae tissue",
"Marine Porifera tissue" against "Sponge tissue", and so on.

The rule now is: **the OMDB value is the label.** Renames are only applied
where the raw string is genuinely ambiguous out of context, and each one is
justified below. Nothing here adds an interpretation the metadata does not
support, which the old dict did:

  "Animal bone biofilm" was displayed as "Whale-fall bone biofilm". All 298
  such genomes come from one study (SUMM22-1) and OMDB places every one of
  them in a "Coastal intertidal zone", not the deep sea where whale falls
  occur. The whale-fall reading was an inference, it was wrong, and it was
  being asserted in a figure. The raw label stands.

  "Seawater" was displayed as "Open water (seawater)". OMDB does not
  distinguish open from coastal water in this field, so "open" was likewise
  an addition. Kept as "Seawater".

EXCLUDE holds the `ecosystem_compartment` values that are not marine habitats
-- freshwater systems, lab microcosms, controls and blanks -- which OMDB files
under the ocean project regardless. It lived in two scripts in two slightly
different versions; one copy now.
"""

# ecosystem_compartment values that are not an ocean habitat. OMDB labels the
# whole database "ocean" at project level, including freshwater and lab work.
EXCLUDE = {
    'NA', '', 'Control', 'Synthetic', 'Freshwater river water', 'Freshwater lake water',
    'Freshwater lake sediment', 'Freshwater aquaculture water', 'Freshwater pond water',
    'Glacier-fed stream water', 'Aquifer water', 'Polar desert soil', 'Seawater microcosm',
    'Sediment microcosm', 'Negative control',
}

# Habitats smaller than this give percentages too unstable to plot or test.
MIN_N = 200

# Only entries whose raw string is unclear standing alone in a figure. Each adds
# information from the taxon name itself, never from an assumption about the site.
RENAME = {
    'Marine Porifera tissue': 'Marine sponge (Porifera) tissue',
    'Marine Hydrozoa tissue': 'Marine hydrozoan tissue',
    'Marine algae thallus': 'Marine algal thallus',
}


def label_for(compartment):
    """Display label for one OMDB ecosystem_compartment value."""
    return RENAME.get(compartment, compartment)
