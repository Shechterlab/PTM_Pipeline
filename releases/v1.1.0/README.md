# Human methylarginine sites, v1.1.0

13,358 arginine methylation sites on 4,722 reviewed human proteins, each with its sequence, methylation evidence, structural context and the studies that report it. This is the site union behind DeAngelo, Silverstein & Shechter, Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, corrected: v1.1.0 replaces v1.0.0 (the published union); CHANGES.md says what changed and why.

## Files

| File | Holds |
|---|---|
| human_methylarginine_sites_v1.1.0.tsv | the release: one row per site |
| human_methylarginine_sites_v1.1.0.xlsx | the same table as a spreadsheet |
| CHANGES.md | what changed since the published version and why |
| supporting/source_records_v1.1.0.tsv | one row per source record: what each study reported and where it was placed |
| supporting/site_changes_v1.0.0_to_v1.1.0.tsv | one row per site of either version: kept, removed or added, and why |
| MD5SUMS | checksums (`md5sum -c MD5SUMS`) |

## Columns

**Site**

- `site_id`: UniProt accession and arginine position, e.g. P35637:R216.
- `gene`: primary gene name of the UniProt entry.
- `protein_name`: recommended protein name of the UniProt entry.
- `uniprot_accession`: reviewed canonical UniProtKB accession (linked to its UniProt page in the workbook).
- `uniprot_entry_name`: UniProtKB entry name.
- `position`: position of the arginine in the reference sequence (1-based).
- `residue`: residue at that position in the reference sequence (always R).

**Sequence**

- `window_15`: residues -7 to +7 around the site; _ beyond a protein end.
- `window_31`: residues -15 to +15 around the site; _ beyond a protein end.
- `protein_length`: length of the reference sequence.
- `reference_sequence_md5`: MD5 of the reference sequence (UniProt sequence MD5 format: uppercase hex); compare with a later UniProt entry to see whether the numbering still holds.

**Methylation evidence**

- `methyl_states_reported`: methylation degree given by the source records: me1 (monomethyl), me2 (dimethyl, asymmetric or symmetric not stated); blank where the sources give none.
- `maron2021_ip_states`: degree called in Maron 2021 Table S2 (PTMScan methyl-arginine immunoprecipitation): me1, me2.
- `maron2021_ip_antibodies`: antibodies whose immunoprecipitation recovered the site in Maron 2021 Table S3: ADMA (asymmetric dimethyl), SDMA (symmetric dimethyl), MMA (monomethyl).
- `maron2021_GSK591_PRMT5i`: response to the type II (PRMT5) inhibitor GSK591 in Maron 2021 Table S3: decreased or increased (adjusted p < 0.05 in at least one immunoprecipitation), mixed, or no significant change.
- `maron2021_MS023_typeIi`: response to the type I PRMT inhibitor MS023 in Maron 2021 Table S3, as for GSK591.
- `maron2021_inhibitor_detail`: per inhibitor and immunoprecipitation: log2 intensity difference vs DMSO and adjusted p value of the strongest peptide (lowest adjusted p). Measured on peptides: arginines on the same peptide share its values.

**Context**

- `motifs`: sequence motifs around the site, as defined by scripts/08: RG, RGG, GRG (G at -1 and +1), GAR (G at -2, A at -1), RXR (R at -2 or +2), CARM1_proline_rich (P within -5..+5), CARM1_hydrophobic_proline (P within -5..+5 and a hydrophobic residue within -3..+3), DR, ER (D or E at -1), RD, RE (D or E at +1).
- `disorder`: MobiDB predicted disorder at the site (scripts/06): disordered; disorder boundary (outside, within 20 residues of a disordered region); ordered; no disorder annotation.
- `disorder_edge_distance`: residues to the nearest edge of a predicted disordered region.
- `domain`: InterPro 108.0 domain or repeat context (scripts/05): in domain; boundary (within 20 residues of one); inter-domain linker; distal; no domain annotation (InterPro has no domain or repeat on the protein).
- `nearest_domain`: InterPro entry nearest the site (the innermost when inside several), with its accession and span.
- `domain_edge_distance`: residues to the nearest domain edge (0 inside a domain).

**Provenance**

- `n_studies`: number of studies or databases reporting the site.
- `studies`: the studies or databases reporting the site; PhosphoSitePlus is listed as PSP, as its terms require (see Sources).
- `study_dois`: their DOIs, in the same order.
- `n_records`: number of source records placed on the site.
- `record_ids`: the records behind the site: record_id in supporting/source_records.

**Status**

- `flags`: cautions, separated by ";" (see Flags).
- `same_window_in`: other reviewed entries holding the same 31 residues (see the shared_sequence flag).

**Change**

- `in_previous_version`: whether the previous version (the published union) had this site.
- `change_since_previous_version`: kept or added; removed sites are listed in supporting/site_changes.
- `change_reasons`: for an added site, why (reason codes, see CHANGES.md).

## Flags

- `pending_splice_junction` (10 sites): a window crossing a splice junction: one full flank (at least 6 residues) matches the entry exactly at the reported arginine, the other does not, and no other entry holds the window. The pipeline leaves it unplaced; the release includes it until the policy is decided.
- `paralog_window` (22 sites): a supporting window matches a paralog exactly and this entry within 2 residues (issue 9, open; accepted in config/remap_evidence_allowlist.tsv).
- `near_match_only` (12 sites): every supporting record matches the entry within 2 residues, not exactly (a variant or sequence conflict).
- `shared_sequence` (593 sites): residues -15 to +15 occur identically in another reviewed entry (listed in same_window_in): histone variants, identical paralogs. The evidence cannot tell the proteins apart; the gene is the one the source named.
- `single_study` (6,452 sites): one study or database reports the site.

## Sources

Sites come from three source tables: Maron et al. 2021 Table S5 (a compilation of ten studies), ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of an iPTMnet-format PTM master table (UniProt and IEDB records). Cite the studies of the sites you use:

| Study or database | Citation | DOI | Records | Sites | Sites only here |
|---|---|---|---|---|---|
| Larsen 2016 | Larsen et al., Sci Signal 9 (2016) | [10.1126/scisignal.aaf7329](https://doi.org/10.1126/scisignal.aaf7329) | 8,881 | 8,500 | 2,618 |
| PhosphoSitePlus (PSP) | Hornbeck et al., Nucleic Acids Res 43:D512-D520 (2015) | [10.1093/nar/gku1267](https://doi.org/10.1093/nar/gku1267) | 8,375 | 8,171 | 2,086 |
| Maron 2021 | Maron et al., iScience 24:102971 (2021) | [10.1016/j.isci.2021.102971](https://doi.org/10.1016/j.isci.2021.102971) | 1,829 | 1,821 | 599 |
| Guo 2014 | Guo et al., Mol Cell Proteomics 13:372-387 (2014) | [10.1074/mcp.O113.027870](https://doi.org/10.1074/mcp.O113.027870) | 1,471 | 1,444 | 68 |
| Musiani 2019 | Musiani et al., Sci Signal 12:eaat8388 (2019) | [10.1126/scisignal.aat8388](https://doi.org/10.1126/scisignal.aat8388) | 1,712 | 1,395 | 253 |
| Fong 2019 | Fong et al., Cancer Cell 36:194-209 (2019) | [10.1016/j.ccell.2019.07.003](https://doi.org/10.1016/j.ccell.2019.07.003) | 1,298 | 1,297 | 182 |
| ProMetheusDB | Massignani et al., Mol Cell Proteomics 21:100243 (2022) | [10.1016/j.mcpro.2022.100243](https://doi.org/10.1016/j.mcpro.2022.100243) | 1,598 | 1,275 | 143 |
| Li 2021 | Li et al., Nat Commun 12:1946 (2021) | [10.1038/s41467-021-21963-1](https://doi.org/10.1038/s41467-021-21963-1) | 1,059 | 1,055 | 180 |
| Fedoriw 2019 | Fedoriw et al., Cancer Cell 36:100-114 (2019) | [10.1016/j.ccell.2019.05.014](https://doi.org/10.1016/j.ccell.2019.05.014) | 951 | 925 | 86 |
| UniProt | UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025) | [10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010) | 762 | 739 | 33 |
| Lim 2020 | Lim et al., Proteome Sci 18:6 (2020) | [10.1186/s12953-020-00162-8](https://doi.org/10.1186/s12953-020-00162-8) | 650 | 602 | 113 |
| Wei 2020 | Wei et al., Front Cell Dev Biol 8:545934 (2020) | [10.3389/fcell.2020.545934](https://doi.org/10.3389/fcell.2020.545934) | 289 | 279 | 55 |
| IEDB | Vita et al., Nucleic Acids Res 47:D339-D343 (2019) | [10.1093/nar/gky1006](https://doi.org/10.1093/nar/gky1006) | 85 | 84 | 36 |

Maron et al. 2021 Tables S2 and S3 annotate sites but add none: methylation degree, the antibodies that recovered each site and its response to PRMT inhibitors (957 sites). A Table S2 or S3 site is attached only when its peptide lies on the reference entry at its stated start, or occurs there once (2,880 of 2,884 Table S2 rows, 24,680 of 24,696 Table S3 site rows).

## How it was built

`bash hpc/build_release.sh v1.1.0` at PTM_Pipeline commit dc2b144940a37c7101654005a952f0ce5e11325d: it parses the source tables, places each record on the reference sequence (scripts/02_remap_arg_methyl_to_canonical.py), stops if any record is placed against its own evidence (scripts/qc_remap_evidence_check.py, 22 accepted paralog rows aside) or any site is not an arginine, then writes this folder (scripts/build_release.py).

- Inputs: inputs/, checked against inputs/MANIFEST.tsv before the build.
- Reference sequences: UniProtKB reviewed human canonical (proteome UP000005640), staged 2026-04-07; the UniProt release was not recorded and predates 2026_03. Positions refer to these sequences; reference_sequence_md5 tells whether a later UniProt entry still has the same sequence. Accession history: UniProt 2026_03.
- Context: MobiDB predicted disorder (export of 2026-03-30), InterPro 108.0 domains and repeats, sequence motifs; the classifiers of scripts/06, scripts/05 and scripts/08.

## Licences

The release is CC BY-NC-SA 4.0 (LICENSE.md). It includes PhosphoSitePlus® (PSP, www.phosphosite.org) records compiled in Maron et al. 2021 Table S5. PhosphoSite, created by Cell Signaling Technology, was licensed under CC BY-NC-SA 3.0 Unported in 2021, when those records were downloaded (inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf). Those terms ask that PSP be cited as Hornbeck et al., Nucleic Acids Res 43:D512-D520 (2015), and that bulk downloads list PSP as the source (the studies column does), and that web pages showing PSP-derived sites link them to the PSP protein page (www.phosphosite.org/uniprotAccAction?id=<UniProt accession>).

## Run record

```
built: 2026-10-20T02:10:32Z
commit: dc2b144940a37c7101654005a952f0ce5e11325d
uncommitted changes:
python: 3.11.17 pandas 2.2.3 numpy 2.2.6 openpyxl 3.1.5
inputs: inputs (all match MANIFEST.tsv)
```
