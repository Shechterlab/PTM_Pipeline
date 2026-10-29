# Changes from v1.4.1 to v1.4.2

| | v1.4.1 | v1.4.2 |
|---|---:|---:|
| Sites | 11,632 | 11,632 |
| Proteins | 4,439 | 4,439 |
| Sites kept pending | 0 | 0 |
| Source records | 67,516 | 67,516 |
| Source records placed | 38,325 | 38,325 |
| Studies and databases | 19 | 19 |

Of the 11,632 sites of v1.4.1 (pending ones included), 11,632 are kept and 0 removed; 0 sites were added. The records of v1.4.1 are read in the same order and all are placed as before. 

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## New annotation layers

v1.4.2 adds a layer of annotation beside the sites. It adds, removes or moves no site; its columns, tables and tests are in `README.md`, `data_dictionary.tsv` and the statistics table, and scripts/release_layers.py builds them.

### Maps of the cell

Where the methylated proteins sit in the cell, by four maps beside the Human Protein Atlas: the Human Cell Map and its finer domains, organelle IP, the nuclear speckle, lamina and PML-body proteomes, and chromatin probability (protein columns `hcm_compartment`, `hcm_domain`, `organelle_ip_compartment`, `nuclear_body_proteomes`, `chromatin_probability`). Each map tests the methylated proteins against the proteins it places, with enricher and within PaxDb abundance deciles. And where the RNAs that proteins with eCLIP data bind gather, by two RNA maps (APEX-seq, nuclear-body APEX2-seq) and two DNA maps (TSA-seq, Hi-C subcompartments), each protein's bound RNAs against the other RNAs a map places within expression and gene-length deciles (`supporting/eclip_rna_compartments`, column `eclip_rna_compartments`); the methylated proteins are compared with the other eCLIP'd proteins in each compartment. The website showed these maps as a preview before this release.

- Compartments where the methylated proteins are enriched at q <= 0.05 among the proteins each map places: Human Cell Map compartments 6 of 20 (4 within abundance deciles); Human Cell Map domains 4 of 23 (4 within abundance deciles); Organelle IP 7 of 20 (7 within abundance deciles); Speckle, lamina and PML-body proteomes 1 of 3 (1 within abundance deciles); Chromatin probability 1 of 4 (1 within abundance deciles).
- RNAs bound by proteins with eCLIP data: of the up to 220 such proteins a map tests, 179 carry a methylarginine. No compartments separate them from the others within abundance deciles (lowest q 0.11), and none by the size of their odds ratios (lowest q 0.23).

Source: Go et al., Nature 595:120 (2021), doi:10.1038/s41586-021-03592-2; Hein et al., Cell 188:1137 (2025), doi:10.1016/j.cell.2024.11.028; Kustatscher et al., EMBO J 33:648 (2014), doi:10.1002/embj.201387614; Barutcu et al., Mol Cell 82:1035 (2022), doi:10.1016/j.molcel.2021.12.010; Fazal et al., Cell 178:473 (2019), doi:10.1016/j.cell.2019.05.027; Chen et al., J Cell Biol 217:4025 (2018), doi:10.1083/jcb.201807108; Rao et al., Cell 159:1665 (2014), doi:10.1016/j.cell.2014.11.021; ENCODE Project Consortium, Nature 583:699 (2020), doi:10.1038/s41586-020-2493-4. Licence: Human Cell Map CC BY 4.0; organelle IP BSD-3; Kustatscher et al. 2014 CC BY; Barutcu et al. 2022 Table S1 (Elsevier): the calls derived from it only; GEO and ENCODE data without restrictions; GENCODE open.

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (Site specificity in the CHANGES.md of v1.2.0): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.4.1_to_v1.4.2.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.4.2.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
