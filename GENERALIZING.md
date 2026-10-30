# Generalising the pipeline and explorer beyond methylarginine

The aim is to extract a codebase that builds the same release and website for another modification (phosphosites, lysine methylation, acetylation …). This note records what is generic already, what ties the code to methylarginine, and an order of work that keeps the site working at every step. New code should follow it now.

## Generic already

- **Release layers.** `scripts/release_layers.py` registers each source as a `Layer` (folder, version it arrived in, site columns, codes, supporting tables, website files, CHANGES text). `ADDING_A_SOURCE.md` describes adding one.
- **Statistics.**
  - `scripts/common.py`: Benjamini–Hochberg, and the Mantel–Haenszel odds ratio with its Robins–Breslow–Greenland interval.
  - `scripts/enrichment.py`: clusterProfiler `enricher` jobs.
  - `scripts/abundance_strata.py`: the PaxDb-decile sensitivity tests.
- **Site-level layers that only need (protein, position).** `scripts/cell_maps.py` and `scripts/interaction_sources.py` work on any site, apart from the constants named below.
- **Website pieces that only need rows.**
  - The chart kit: `hbars`, `chartSvg`, `chartTsv`, `figExportSvg`, `withCredit`, `exportFigure`.
  - The tile kit: `aHome`, `aForest`, `aMeters`.
  - The tests on a selection: `selSetTest`, `mhTest`, `hyperSf`, `bhAdjust`.
  - How each test is described in exports: `TEST_WORDS`, with the counted nouns in `WORDS`.

## What ties it to methylarginine (counted 30 October 2026)

- **The explorer is one file.** `docs/assets/explorer.js` has 17,600 lines in one closure. "methylarginine" is on 461 of its lines, "arginine" on 699 and "PRMT" on 108.
- **The residue is a literal.** 27 pipeline scripts compare residues with `'R'` directly (for example `02_remap_arg_methyl_to_canonical.py`, `08_motif_and_arg_odds.py`, `10_methyl_arg_clustering.py` and `37_export_interaction_data.py`).
- **Mark-specific vocabulary:**
  - methylation states (`Rme1`, `Rme2a`, `Rme2s`);
  - PRMT types and what they make (`PRMT_TYPE` and `PRMT_MAKES` in `interaction_sources.py`);
  - Tudor readers;
  - the PSI-MOD terms and PDB chemical components of a methylated arginine (`METHYL_R`: 2MR, DA2, NMM);
  - the motif library (RG, RGG, GAR, RXR, the CARM1 proline motifs …);
  - the PRMT inhibitor and knockout studies.
- **Names and texts.**
  - Data files are named `methylarg_*.json`.
  - The page title, version texts and citation.

## Proposal

1. **A modification profile.** One JSON file per modification, read by the pipeline and the website. It holds:
   - the residue letter(s);
   - the mark's states and their codes;
   - the enzyme families (writers, erasers, readers), as UniProt lists;
   - the PSI-MOD terms and PDB chemical components;
   - the motifs;
   - the nouns used in text (`WORDS`);
   - the site's name and the paper to cite.
2. **A pipeline package.** Staging helpers, remapping to canonical sequences, site tests, the layer registry, enrichment, abundance strata, cell maps and interaction provenance move into one Python package, called with a profile. Sources specific to one modification stay as layers (the Maron parsers, the PRMT perturbation studies).
3. **An explorer split into ES modules, along seams it already has.** Behaviour stays the same; no framework and no build step are needed beyond, at most, one bundling command.
   - formatting and utilities;
   - statistics (`mhTest`, `hyperSf`, `fisherP`, `mwuP`, `bhAdjust`);
   - exports (`chartSvg`, `figExportSvg`, `withCredit`, `exportFigure`, `TEST_WORDS`);
   - the chart and tile kits;
   - state, URL and history;
   - the pages (Overview, Proteins, Sites, Interactions, protein page);
   - panels specific to one modification, registered the way release layers are (perturbations, writers and readers).
4. **A written data contract.** JSON schemas for the files in `docs/data/`, versioned with the release. File names drop the `methylarg_` prefix, through a mapping that keeps old releases readable.

## Order of work

Each step leaves the site working, checked by the same browser tests.

1. **Freeze the contract.** Write the schemas for `docs/data/` and a Playwright smoke test of the main paths (Overview, a selection, a protein page, a site card, the Interactions page, Back and Forward, every export button).
2. **Words first.** Route the nouns through `WORDS` and the profile, file by file. This is begun: the test wording of every exported figure comes from `TEST_WORDS`.
3. **Split `explorer.js` into modules** along the seams above, one at a time, with the smoke test after each.
4. **Make the pipeline profile-driven.** Replace the `'R'` literals with the profile's residue, and move the shared code into the package.
5. **Pilot a second modification.** Lysine methylation is a good candidate, since CPLM4 is already staged by `scripts/39_stage_cplm4_lys_methylation.py`. The pilot will find the couplings this list misses.
