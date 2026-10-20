# Human methylarginine sites, v1.1.0

13,358 sites on 4,722 proteins, and 10 pending sites, 13,368 rows: the human protein arginine methylation sites of DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with corrections to how some source records were placed on UniProt sequences. [CHANGES.md](CHANGES.md) says what changed since v1.0.0 (the published version) and why.

Every site has its sequence, methylation evidence, structural context, the studies that report it, and the source records behind it. Positions are on UniProt release 2026_03 (reviewed human reference proteome UP000005640).

## Files

| File | One row per |
|---|---|
| `human_methylarginine_sites_v1.1.0.tsv` | site, including the 10 sites kept pending (flag `one_sided_window`) |
| `human_methylarginine_sites_v1.1.0.xlsx` | the same table as a spreadsheet |
| `human_methylarginine_proteins_v1.1.0.tsv` | protein: names, site counts, functional classes, arginine enrichment, site density and clusters, condensates, locations, disease families, other PTMs, interactions |
| `studies_v1.1.0.tsv` | study or database: citation, DOI, terms, counts |
| `statistics_v1.1.0.tsv` | enrichment test behind the summaries: classes, motifs, disorder, domains, disease families, condensates, locations, interactions |
| `supporting/source_records_v1.1.0.tsv` | source record: what the study reported, its evidence, where each version placed it |
| `supporting/site_changes_v1.0.0_to_v1.1.0.tsv` | site of either version: kept, removed or added, and why |
| `supporting/maron2021_detection_v1.1.0.tsv` | site detected in Maron et al. 2021 Table S2 or S3 (per antibody) |
| `supporting/maron2021_inhibitors_v1.1.0.tsv` | site, antibody and PRMT inhibitor in Maron et al. 2021 Table S3 |
| `supporting/carm1_ko_v1.1.0.tsv` | site and cell line in the CARM1 knockout data of Shishkova et al. 2017 |
| `supporting/site_interactions_v1.1.0.tsv` | piece of interaction evidence at a site: crosslinks, PDB interfaces, predicted interfaces and contacts, IntAct records |
| `supporting/residue_evidence_v1.1.0.tsv` | run of residues of a protein with the same interaction evidence |
| `supporting/protein_partners_v1.1.0.tsv` | protein and STRING physical partner (score >= 0.7) |
| `supporting/complexes_v1.1.0.tsv` | Complex Portal or hu.MAP complex with a methylated member |
| `supporting/disease_terms_v1.1.0.tsv` | enriched DisGeNET term (DOSE enrichDGN) |
| `supporting/disease_families_v1.1.0.tsv` | disease family: terms and genes |
| `supporting/density_windows_v1.1.0.tsv` | protein and window width with significantly dense sites (q <= 0.05) |
| `supporting/condensates_v1.1.0.tsv` | CD-CODE condensate: name, type, members and methylated members |
| `supporting/hpa_locations_v1.1.0.tsv` | protein and Human Protein Atlas location |
| `supporting/other_ptm_sites_v1.1.0.tsv` | other modification on a protein of the release (PTM master table) |
| `data_dictionary.tsv` | column of every table, with its codes and their meanings |
| `release.json` | the release: version, counts, reference, inputs, files with MD5 |
| `CHANGES.md` | what changed since the previous version |
| `MD5SUMS` | checksums (`md5sum -c MD5SUMS`) |

Tables are UTF-8, tab-separated, one header row. Lists are separated by `; `, booleans are `true`/`false`, and a missing value is empty. Coded columns hold codes only; `data_dictionary.tsv` gives their meanings.

## Cite

Cite DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, and the studies of the sites you use (`studies` column, `studies_v1.1.0.tsv`):

| Study or database | Citation | DOI | Records | Sites | Only here |
|---|---|---|---:|---:|---:|
| Larsen 2016 (`Larsen2016`) | Larsen SC et al., Sci Signal 9:rs9 (2016) | [10.1126/scisignal.aaf7329](https://doi.org/10.1126/scisignal.aaf7329) | 8,881 | 8,491 | 2,609 |
| PhosphoSitePlus (PSP) | Hornbeck PV et al., Nucleic Acids Res 43:D512-D520 (2015) | [10.1093/nar/gku1267](https://doi.org/10.1093/nar/gku1267) | 8,375 | 8,171 | 2,086 |
| Maron 2021 (`Maron2021`) | Maron MI et al., iScience 24:102971 (2021) | [10.1016/j.isci.2021.102971](https://doi.org/10.1016/j.isci.2021.102971) | 1,829 | 1,821 | 599 |
| Guo 2014 (`Guo2014`) | Guo A et al., Mol Cell Proteomics 13:372-387 (2014) | [10.1074/mcp.O113.027870](https://doi.org/10.1074/mcp.O113.027870) | 1,471 | 1,443 | 67 |
| Musiani 2019 (`Musiani2019`) | Musiani D et al., Sci Signal 12:eaat8388 (2019) | [10.1126/scisignal.aat8388](https://doi.org/10.1126/scisignal.aat8388) | 1,712 | 1,395 | 253 |
| Fong 2019 (`Fong2019`) | Fong JY et al., Cancer Cell 36:194-209 (2019) | [10.1016/j.ccell.2019.07.003](https://doi.org/10.1016/j.ccell.2019.07.003) | 1,298 | 1,297 | 182 |
| ProMetheusDB | Massignani E et al., Mol Cell Proteomics 21:100243 (2022) | [10.1016/j.mcpro.2022.100243](https://doi.org/10.1016/j.mcpro.2022.100243) | 1,598 | 1,275 | 143 |
| Li 2021 (`Li2021`) | Li WJ et al., Nat Commun 12:1946 (2021) | [10.1038/s41467-021-21963-1](https://doi.org/10.1038/s41467-021-21963-1) | 1,059 | 1,055 | 180 |
| Fedoriw 2019 (`Fedoriw2019`) | Fedoriw A et al., Cancer Cell 36:100-114 (2019) | [10.1016/j.ccell.2019.05.014](https://doi.org/10.1016/j.ccell.2019.05.014) | 951 | 925 | 86 |
| UniProt | UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025) | [10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010) | 762 | 739 | 33 |
| Lim 2020 (`Lim2020`) | Lim Y et al., Proteome Sci 18:6 (2020) | [10.1186/s12953-020-00162-8](https://doi.org/10.1186/s12953-020-00162-8) | 650 | 602 | 113 |
| Wei 2020 (`Wei2020`) | Wei M et al., Front Cell Dev Biol 8:545934 (2020) | [10.3389/fcell.2020.545934](https://doi.org/10.3389/fcell.2020.545934) | 289 | 279 | 55 |
| IEDB | Vita R et al., Nucleic Acids Res 47:D339-D343 (2019) | [10.1093/nar/gky1006](https://doi.org/10.1093/nar/gky1006) | 85 | 84 | 36 |

The site records come from three source tables: Maron et al. 2021 Table S5 (sites compiled from ten studies), ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of an iPTMnet-format PTM master table (UniProt and IEDB records). Maron et al. 2021 Tables S2 and S3 and the CARM1 knockout data (Shishkova et al. 2017, doi:10.1038/ncomms15571) annotate sites but add none.

## Licence

The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)); each study keeps its own terms (`licence` in the studies file). The release includes PhosphoSitePlus® (PSP, www.phosphosite.org) records compiled in Maron et al. 2021 Table S5. PhosphoSite, created by Cell Signaling Technology, was licensed under CC BY-NC-SA 3.0 Unported in 2021, when those records were downloaded ([inputs/licences/](../../inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf)). Those terms ask that PSP be cited (Hornbeck et al., Nucleic Acids Res 43:D512-D520, 2015), that bulk downloads list PSP as the source (the `studies` column does), and that web pages showing PSP-derived sites link them to the PSP protein page (www.phosphosite.org/uniprotAccAction?id=<UniProt accession>).

The annotation and interaction sources keep their own terms too (`licence` in `inputs/MANIFEST.tsv`). Three are non-commercial: the DisGeNET terms behind the disease tables (CC BY-NC-SA 4.0), CD-CODE condensates (CC BY-NC 4.0) and the pepR-MS RNA crosslink compendium (CC BY-NC).

## How it was built

`bash hpc/build_release.sh v1.1.0` at commit `9af7bb9550f5a73dc9af85f16a6777a268e59163` of https://github.com/Shechterlab/PTM_Pipeline:

1. checks every input against `inputs/MANIFEST.tsv`;
2. parses the source tables and places each record on the reference sequence (`scripts/02_remap_arg_methyl_to_canonical.py`);
3. stops if any record is placed against its own evidence (`scripts/qc_remap_evidence_check.py`; the 22 accepted paralog rows of `config/remap_evidence_allowlist.tsv` aside) or any site is not an arginine;
4. runs the analyses on the union: functional classes (`scripts/03`), motifs and arginine enrichment (`08`), disorder and domain context (`06`, `05`, `07`), other PTMs (`09`), site density and clusters (`10`, `28`, `25`), condensates (`12`), Maron et al. 2021 summaries (`16`) and disease enrichment (`14`);
5. encodes the release for the website (`scripts/36`) and computes the interaction evidence of its sites from `inputs/interactions/` (`scripts/37`);
6. writes this folder (`scripts/build_release.py`).

- **Reference sequences:** UniProt release 2026_03 (02-September-2026), the 20,416 reviewed entries of the human reference proteome UP000005640 (`inputs/reference/`). `reference_sequence_md5` tells whether a later UniProt entry still has the same numbering.
- **Context:** MobiDB predicted disorder (export of 2026-03-30), InterPro 108.0 domains and repeats, sequence motifs: the classifiers of `scripts/06`, `scripts/05`, `scripts/07` and `scripts/08`. The intervals were computed on earlier UniProt sequences and moved to these (`inputs/annotation/interval_lift_to_uniprot_2026_03.tsv`).
- **Previous version:** `releases/v1.0.0/` (the published union, git tag v1.0.0).

## Columns

The same definitions, with every code and its meaning, are in `data_dictionary.tsv`.

### `human_methylarginine_sites_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | UniProt accession and arginine position, e.g. P35637:R216. |
| `uniprot_accession` | string | reviewed UniProtKB accession of the canonical sequence. |
| `uniprot_entry_name` | string | UniProtKB entry name. |
| `gene` | list&lt;string&gt; | primary gene name of the UniProt entry (several for a few entries, e.g. histone genes), empty when UniProt gives none. |
| `protein_name` | string | recommended protein name of the UniProt entry. |
| `position` | integer | position of the arginine in the reference sequence (1-based, initiator methionine = 1). |
| `residue` | code | residue at that position in the reference sequence: always R. |
| `window_15` | string | residues -7 to +7 around the site, _ beyond a protein end. |
| `window_31` | string | residues -15 to +15 around the site, _ beyond a protein end. |
| `protein_length` | integer | length of the reference sequence. |
| `reference_sequence_md5` | string | MD5 of the reference sequence, uppercase hex as UniProt gives it (a later UniProt entry with the same MD5 has the same numbering). |
| `methyl_states_reported` | list&lt;code&gt; | methylation degrees given by the source records, in code order, empty where they give none. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`. |
| `maron2021_s2_states` | list&lt;code&gt; | degrees called at the site in Maron et al. 2021 Table S2. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`. |
| `maron2021_s3_antibodies` | list&lt;code&gt; | antibodies whose immunoprecipitation recovered the site in Maron et al. 2021 Table S3 (supporting/maron2021_detection). Codes: `MMA`, `ADMA`, `SDMA`. |
| `maron2021_GSK591_call` | code | response to GSK591 in Maron et al. 2021 Table S3, combined over antibodies (supporting/maron2021_inhibitors: a peptide changes at adjusted p <= 0.05), empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_change`. |
| `maron2021_MS023_call` | code | response to MS023, as for GSK591. Codes: `decreased`, `increased`, `mixed`, `no_change`. |
| `carm1_ko_call` | code | response to CARM1 knockout (Shishkova et al. 2017), combined over cell lines (supporting/carm1_ko), empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_change`. |
| `motifs` | list&lt;code&gt; | sequence motifs around the site (scripts/08), in code order. Codes: `RG`, `RGG`, `GRG`, `GAR`, `RXR`, `CARM1_proline_rich`, `CARM1_hydrophobic_proline`, `DR`, `RD`, `ER`, `RE`. |
| `disorder` | code | MobiDB predicted disorder at the site (scripts/06). Codes: `disordered`, `disorder_boundary`, `ordered`, `no_disorder_annotation`. |
| `disorder_edge_distance` | integer | residues to the nearest edge of a disordered region, empty with no annotation. |
| `domain` | code | InterPro 108.0 domain or repeat context (scripts/05). Codes: `in_domain`, `boundary`, `inter_domain_linker`, `distal`, `no_domain_annotation`. |
| `nearest_domain_interpro` | string | InterPro accession of the nearest domain or repeat (the innermost when inside several). |
| `nearest_domain_name` | string | its InterPro name. |
| `nearest_domain_start` | integer | its first residue. |
| `nearest_domain_end` | integer | its last residue. |
| `nearest_domain_class` | string | its class in config/domain_class_ontology.json, as scripts/07 assigns it, No domain annotation when none. |
| `domain_side` | code | where the site lies relative to the nearest domain, empty when none. Codes: `inside`, `N_terminal_to_domain`, `C_terminal_to_domain`. |
| `domain_edge_distance` | integer | residues to the nearest domain edge (0 inside a domain), empty when none. |
| `interaction_status` | code | interaction evidence at the site (supporting/site_interactions), empty for the sites kept pending. Codes: `experimental`, `predicted_only`, `none`. |
| `interaction_evidence` | list&lt;code&gt; | the evidence types at the site, in code order. Codes: `rna_crosslink`, `rbdmap_fragment`, `dna_crosslink`, `dna_binding_region`, `pdb_protein`, `pdb_rna`, `pdb_dna`, `pdb_self`, `predicted_interface`, `predicted_contact`, `intact_mutation`, `intact_binding_region`, `intact_methylated`. |
| `n_studies` | integer | number of studies or databases reporting the site. |
| `studies` | list&lt;code&gt; | the studies or databases reporting the site, in the order of the studies file, PhosphoSitePlus is PSP, as its terms ask. Codes: `Larsen2016`, `PSP`, `Maron2021`, `Guo2014`, `Musiani2019`, `Fong2019`, `ProMetheusDB`, `Li2021`, `Fedoriw2019`, `UniProt`, `Lim2020`, `Wei2020`, `IEDB`. |
| `study_dois` | list&lt;string&gt; | their DOIs, in the same order. |
| `source_tables` | list&lt;code&gt; | the source tables of its records, in code order. Codes: `maron2021_s5`, `prometheusdb_mmc2`, `ptm_master`. |
| `n_records` | integer | number of source records placed on the site. |
| `record_ids` | list&lt;integer&gt; | the records behind the site (record_id in supporting/source_records), ascending. |
| `reported_positions` | list&lt;integer&gt; | positions the records report, ascending, only when one differs from position. |
| `flags` | list&lt;code&gt; | cautions, in code order. Codes: `one_sided_window`, `paralog_window`, `near_match_only`, `shared_sequence`, `single_study`. |
| `same_window_accessions` | list&lt;string&gt; | every other reviewed entry holding the 31 residues of window_31, sorted. |
| `same_window_genes` | list&lt;string&gt; | the first primary gene of each (the accession when UniProt gives no gene), same order. |
| `in_previous_version` | boolean | whether the previous version had this site. |
| `change_since_previous_version` | code | kept or added (removed sites are in supporting/site_changes). Codes: `kept`, `added`, `removed`. |
| `change_reasons` | list&lt;code&gt; | why the site was added, or renumbered (CHANGES.md). Codes: `peptide_recentred`, `unreviewed_accession`, `lysine_window`, `truncated_window`, `multi_arginine_peptide`, `peptide_start`, `wrong_protein`, `outdated_accession_or_symbol`, `evidence_mismatch`, `window_padding`, `ambiguous_evidence`, `one_sided_window`, `reference_renumbered`, `other`. |

### `human_methylarginine_proteins_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | reviewed UniProtKB accession. |
| `uniprot_entry_name` | string | UniProtKB entry name. |
| `gene` | list&lt;string&gt; | primary gene name of the UniProt entry (several for a few entries). |
| `gene_synonyms` | list&lt;string&gt; | other gene names UniProt gives the entry (release 2026_03). |
| `reported_genes` | list&lt;string&gt; | other gene symbols the source records used, most frequent first. |
| `protein_name` | string | recommended protein name. |
| `protein_length` | integer | length of the reference sequence. |
| `reference_sequence_md5` | string | MD5 of the reference sequence, uppercase hex as UniProt gives it. |
| `n_sites` | integer | methylarginine sites, not counting the sites kept pending. |
| `n_pending` | integer | sites kept pending (flag one_sided_window). |
| `studies` | list&lt;code&gt; | studies or databases reporting its sites, in the order of the studies file. Codes: `Larsen2016`, `PSP`, `Maron2021`, `Guo2014`, `Musiani2019`, `Fong2019`, `ProMetheusDB`, `Li2021`, `Fedoriw2019`, `UniProt`, `Lim2020`, `Wei2020`, `IEDB`. |
| `functional_classes` | list&lt;string&gt; | broad functional classes (config/functional_ontology.json, scripts/03). |
| `functional_subclasses` | list&lt;string&gt; | functional subclasses. |
| `arginine_count` | integer | arginines in the reference sequence. |
| `arg_enrichment_log2` | number | shrinkage estimate of the log2 enrichment of methylarginines over the proteome-wide rate (scripts/08). |
| `arg_enrichment_q` | number | q value of the per-protein test (scripts/08), empty when not tested. |
| `nearest_spacing_min` | number | smallest distance in residues from a site to the nearest other site (scripts/10). |
| `nearest_spacing_median` | number | median of these distances. |
| `max_sites_3aa` | integer | most sites in any 3-residue window (scripts/10). |
| `max_sites_5aa` | integer | most sites in any 5-residue window. |
| `max_sites_10aa` | integer | most sites in any 10-residue window. |
| `max_sites_20aa` | integer | most sites in any 20-residue window. |
| `max_sites_50aa` | integer | most sites in any 50-residue window. |
| `best_window_aa` | integer | width of the protein's most significant density window: scripts/28 compares the most sites in any window of 5 to 100 residues with a null that places the sites on the protein's own arginines at random (1,000 permutations, proteins with at least 3 sites), and this is the window with the lowest q value (then the largest excess). |
| `best_window_start` | integer | its first residue. |
| `best_window_end` | integer | its last residue. |
| `best_window_sites` | integer | sites in it. |
| `best_window_log2_obs_null` | number | log2 of observed over null mean sites. |
| `best_window_q` | number | its q value. |
| `n_dense_windows` | integer | windows with q <= 0.05 (supporting/density_windows). |
| `architecture_class` | string | domain-cluster architecture of the clustered disordered sites (scripts/25: at least 3 sites within 20 residues), empty when none. |
| `architecture_clusters` | list&lt;string&gt; | the site positions of each cluster, comma-separated. |
| `condensates` | list&lt;string&gt; | CD-CODE condensates the protein belongs to (CD-CODE ids, export of 2026-03-18, scripts/12). |
| `hpa_locations` | list&lt;string&gt; | Human Protein Atlas 25.1 subcellular locations, main locations first (supporting/hpa_locations). |
| `disease_families` | list&lt;string&gt; | disease families whose genes include the protein's gene (scripts/14). |
| `disease_family_terms` | list&lt;integer&gt; | for each of these families, how many of its DisGeNET terms name the gene. |
| `als_ftd_family_gene` | boolean | the gene is in the ALS / motor neuron / FTD family. |
| `n_phospho` | integer | phosphorylation sites in the PTM master table (supporting/other_ptm_sites). |
| `n_lys_acyl` | integer | lysine acylation sites. |
| `n_ub_sumo` | integer | ubiquitination or SUMOylation sites. |
| `n_lys_methyl` | integer | lysine methylation sites. |
| `n_sites_experimental` | integer | sites with interaction_status experimental. |
| `n_sites_predicted_only` | integer | sites with interaction_status predicted_only. |
| `go_rna_binding` | code | Gene Ontology RNA binding (GO:0003723 or a child term). Codes: `experimental`, `other_evidence`, `none`. |
| `go_dna_binding` | code | Gene Ontology DNA binding (GO:0003677 or a child term). Codes: `experimental`, `other_evidence`, `none`. |
| `n_string_partners` | integer | STRING 12.5 physical partners at score >= 0.7 (supporting/protein_partners). |
| `n_partners_methylated` | integer | of these, proteins of the site table. |
| `n_complexes` | integer | Complex Portal and hu.MAP complexes the protein is in (supporting/complexes). |

### `studies_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `study_id` | code | study or database. Codes: `Larsen2016`, `PSP`, `Maron2021`, `Guo2014`, `Musiani2019`, `Fong2019`, `ProMetheusDB`, `Li2021`, `Fedoriw2019`, `UniProt`, `Lim2020`, `Wei2020`, `IEDB`. |
| `label` | string | name to show. |
| `source_table` | list&lt;code&gt; | the source tables its records come from. Codes: `maron2021_s5`, `prometheusdb_mmc2`, `ptm_master`. |
| `citation` | string | citation, checked against Crossref. |
| `doi` | string | DOI. |
| `licence` | string | terms of the records as this release received them. |
| `n_records` | integer | source records. |
| `n_records_placed` | integer | of these, placed on a site. |
| `n_sites` | integer | sites (status site) it reports. |
| `n_sites_only_here` | integer | sites no other study or database reports. |
| `note` | string |  |

### `statistics_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `analysis` | code | what was tested. Codes: `functional_class`, `functional_subclass`, `motif`, `disorder`, `domain_context`, `domain_class`, `disease_family`, `condensate`, `hpa_location`, `pdb_interface`, `predicted_interface`, `crosslink`, `intact`, `go_binding`. |
| `category` | string | the class, motif, family, condensate (CD-CODE id) or location tested. |
| `stratum` | string | all, or the subset the test is restricted to. |
| `a` | integer | targets in the category. |
| `b` | integer | targets not in the category. |
| `c` | integer | background in the category. |
| `d` | integer | background not in the category. |
| `odds_ratio` | number | odds ratio (method). |
| `log2_odds_ratio` | number | its log2. |
| `p` | number | p value. |
| `q` | number | Benjamini-Hochberg adjusted p value within the analysis. |
| `method` | code | the test. Codes: `fisher`, `mantel_haenszel`. |
| `background` | string | what the targets are compared with. |
| `n_strata` | integer | strata of a stratified test, empty otherwise. |

### `supporting/source_records_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `record_id` | integer | row number of the record in the source tables as the pipeline reads them, the same in every release. |
| `source_table` | code | the source table the record comes from. Codes: `maron2021_s5`, `prometheusdb_mmc2`, `ptm_master`. |
| `study` | code | the study or database that reported it. Codes: `Larsen2016`, `PSP`, `Maron2021`, `Guo2014`, `Musiani2019`, `Fong2019`, `ProMetheusDB`, `Li2021`, `Fedoriw2019`, `UniProt`, `Lim2020`, `Wei2020`, `IEDB`. |
| `study_doi` | string | its DOI. |
| `reference_text` | string | the reference as the source table gives it. |
| `pmid` | string | PubMed ID, where the source gives one. |
| `reported_accession` | list&lt;string&gt; | accession as reported (a few records give two). |
| `reported_gene` | string | gene as reported. |
| `reported_position` | string | position as reported. |
| `reported_methyl_state` | string | methylation state as reported, verbatim. |
| `reported_enzyme` | string | methyltransferase as reported (PTM master table only). |
| `evidence_kind` | code | kind of sequence evidence. Codes: `window_13`, `window_15`, `window_31`, `peptide`, `marked_peptide`, `position_only`. |
| `evidence_sequence` | list&lt;string&gt; | the window, peptide or marked peptide, verbatim (a few records give several). |
| `placed_site` | string | the site this version places the record on, empty when not placed. |
| `placement_code` | code | how the record was placed (scripts/02 remap_status). Codes: see `data_dictionary.tsv`. |
| `evidence_check` | code | independent check of the placement against the record's own evidence (scripts/qc_remap_evidence_check.py). Codes: `agrees`, `agrees_no_site`, `one_of_several`, `ambiguous`, `dropped_evidence_has_site`, `contradicts`. |
| `accepted_exception` | string | why a contradiction is accepted for now (config/remap_evidence_allowlist.tsv). |
| `pending_site` | string | the site a one-sided window record supports (flag one_sided_window), empty otherwise. |
| `previous_site` | string | the site the previous version placed the record on, in its numbering. |
| `change_since_previous_version` | code | how the placement changed. Codes: `same`, `moved`, `newly_placed`, `no_longer_placed`, `not_placed`. |
| `reason_code` | code | why it changed. Codes: `peptide_recentred`, `unreviewed_accession`, `lysine_window`, `truncated_window`, `multi_arginine_peptide`, `peptide_start`, `wrong_protein`, `outdated_accession_or_symbol`, `evidence_mismatch`, `window_padding`, `ambiguous_evidence`, `one_sided_window`, `reference_renumbered`, `other`. |

### `supporting/site_changes_v1.0.0_to_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site: in this version's numbering when the site is in this version. |
| `previous_site_id` | string | the same site in the previous version's numbering, empty for an added site. |
| `gene` | list&lt;string&gt; | primary gene name (for an entry not in the reference, the gene its records report). |
| `uniprot_accession` | string | accession. |
| `position` | integer | arginine position, in this version's numbering when the site is in this version. |
| `in_previous_version` | boolean | a site of the previous version. |
| `in_this_version` | boolean | a site of this version. |
| `change` | code | kept, removed or added. Codes: `kept`, `added`, `removed`. |
| `reason_codes` | list&lt;code&gt; | why it was removed, added or renumbered. Codes: `peptide_recentred`, `unreviewed_accession`, `lysine_window`, `truncated_window`, `multi_arginine_peptide`, `peptide_start`, `wrong_protein`, `outdated_accession_or_symbol`, `evidence_mismatch`, `window_padding`, `ambiguous_evidence`, `one_sided_window`, `reference_renumbered`, `other`. |
| `records_now_at_or_previously_at` | list&lt;string&gt; | for a removed site, where its records are placed now, for an added site, where they were placed before. |
| `record_ids` | list&lt;integer&gt; | the records involved (supporting/source_records). |

### `supporting/maron2021_detection_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `maron_table` | code | the table. Codes: `S2`, `S3`. |
| `antibody` | code | the immunoprecipitation (Table S3 only). Codes: `MMA`, `ADMA`, `SDMA`. |
| `methyl_states` | list&lt;code&gt; | the degrees called (Table S2 only). Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`. |
| `n_peptides` | integer | distinct peptides (sequence and modifications) carrying the site. |

### `supporting/maron2021_inhibitors_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `antibody` | code | the immunoprecipitation. Codes: `MMA`, `ADMA`, `SDMA`. |
| `inhibitor` | code | the treatment, compared with DMSO. Codes: `GSK591`, `MS023`. |
| `n_peptides` | integer | peptide rows carrying the site with a measured difference. |
| `n_decreased` | integer | of these, lower with the inhibitor at adjusted p <= 0.05. |
| `n_increased` | integer | higher with the inhibitor at adjusted p <= 0.05. |
| `log2_diff_median` | number | median log2 intensity difference, inhibitor minus DMSO. |
| `log2_diff_strongest` | number | the difference of the peptide with the lowest adjusted p. |
| `adj_p_min` | number | the lowest adjusted p value (Maron et al. 2021). |
| `call` | code | decreased, increased, mixed or no_change over the peptides, peptides with several arginines give each the same values. Codes: `decreased`, `increased`, `mixed`, `no_change`. |

### `supporting/carm1_ko_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `cell_line` | code | the cell line. Codes: `MCF7`, `MDA-MB-231`. |
| `n_peptides` | integer | peptides carrying the site. |
| `log2_ko_over_wt` | number | median log2 knockout / parental intensity over the peptides. |
| `p` | number | the lowest p value over the peptides (Shishkova et al. 2017). |
| `n_decreased` | integer | peptides at log2 <= -1 and p < 0.05. |
| `n_increased` | integer | peptides at log2 >= 1 and p < 0.05. |
| `call` | code | a peptide decreases at log2 <= -1 and p < 0.05, increases at log2 >= 1 and p < 0.05 (scripts/41), combined over peptides. Codes: `decreased`, `increased`, `mixed`, `no_change`. |

### `supporting/site_interactions_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table). |
| `evidence_type` | code | the kind of evidence. Codes: `rna_crosslink`, `rbdmap_fragment`, `dna_crosslink`, `dna_binding_region`, `pdb_protein`, `pdb_rna`, `pdb_dna`, `pdb_self`, `predicted_interface`, `predicted_contact`, `intact_mutation`, `intact_binding_region`, `intact_methylated`. |
| `localisation` | code | where the evidence falls relative to the site. Codes: `at_residue`, `candidate_residue`, `within_3`, `in_fragment`, `in_peptide`, `in_region`. |
| `partner_kind` | code | what the arginine contacts, empty when the record names no partner. Codes: `protein`, `self`, `rna`, `dna`, `other`. |
| `partner_accession` | string | the partner protein (reviewed UniProtKB accession). |
| `partner_gene` | string | its primary gene name. |
| `source` | string | PDBe-KB, ProtVar, AFDB (AlphaFold DB), IntAct, RBDmap, Trendel2025, UniProt, or the RNA crosslink dataset (Kramer2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS, pepR-MS_focused). |
| `structures` | list&lt;string&gt; | PDB entries (PDB rows), or the predicted model (ProtVar pair or AlphaFold DB model id). |
| `score` | number | pDockQ (ProtVar) or ipTM (AlphaFold DB) of the predicted complex. |
| `distance_a` | number | closest distance in the predicted model between the arginine and the partner, in angstroms. |
| `detail` | string | bonds, pLDDT and partner residues of a predicted contact; effect, residues, PubMed and IntAct ids of an IntAct record. |
| `counts_as_evidence` | boolean | the record is site-level evidence: crosslinks at or as a candidate residue, RBDmap fragments, PDB interfaces, predicted interfaces, confident predicted contacts, and IntAct mutations that change an interaction, short binding regions and methylated residues. interaction_status counts all of these except IntAct. |

### `supporting/residue_evidence_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table. |
| `start` | integer | first residue of a run of consecutive residues with the same evidence. |
| `end` | integer | last residue. |
| `evidence_type` | code | the kind of evidence. Codes: `rna_crosslink`, `rbdmap_fragment`, `dna_crosslink`, `dna_binding_region`, `pdb_protein`, `pdb_rna`, `pdb_dna`, `pdb_self`, `predicted_interface`, `predicted_contact`, `intact_mutation`, `intact_binding_region`, `intact_methylated`. |
| `localisation` | code | at the residues, or a region or peptide that holds them. Codes: `at_residue`, `candidate_residue`, `within_3`, `in_fragment`, `in_peptide`, `in_region`. |
| `n_partners` | integer | partners: most along the run (PDB protein, predicted interface) or records (IntAct). |
| `best_score` | number | best pDockQ along the run (predicted interfaces). |
| `sources` | list&lt;string&gt; | sources of the evidence. |

### `supporting/protein_partners_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table. |
| `partner_accession` | string | STRING 12.5 physical partner at score >= 0.7 (reviewed UniProtKB accession). |
| `partner_gene` | string | its primary gene name. |
| `string_score` | number | STRING physical interaction score (0-1). |
| `partner_methylated` | boolean | the partner is a protein of the site table. |
| `partner_classes` | list&lt;code&gt; | reader or writer classes of the partner. Codes: `prmt`, `tudor`, `wd40`, `phd`. |

### `supporting/complexes_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `complex_id` | string | Complex Portal id. |
| `complex_name` | string | its name (empty for most predicted complexes). |
| `curation` | code | curated or predicted. Codes: `curated`, `predicted`. |
| `n_members` | integer | member slots (alternative members of a slot count once). |
| `n_members_methylated` | integer | slots with a protein of the site table. |
| `methylated_members` | list&lt;string&gt; | those proteins. |
| `fold_enrichment` | number | methylated share of the members over the methylated share of reviewed human proteins. |
| `p` | number | hypergeometric p value, P(X >= n_members_methylated). |
| `q` | number | Benjamini-Hochberg over all complexes of the same curation, those without a methylated member included. Only complexes with a methylated member are listed. |

### `supporting/disease_terms_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `term_id` | string | DisGeNET term (UMLS CUI). |
| `term_name` | string | its name. |
| `family` | string | its disease family (scripts/disease_families.py), empty when none. |
| `n_genes_term` | integer | genes of the term in the universe. |
| `n_methylated` | integer | of these, genes of methylated proteins. |
| `fold_enrichment` | number | observed over expected. |
| `p` | number | hypergeometric p value (DOSE enrichDGN). |
| `p_adjust` | number | Benjamini-Hochberg adjusted p value over the tested terms. |
| `q` | number | q value (qvalue, lambda 0.05, as DOSE). |
| `genes` | list&lt;string&gt; | the methylated genes in the term (symbols, org.Hs.eg.db). |
| `uniprot_accessions` | list&lt;string&gt; | the proteins of the site table those genes stand for. |

### `supporting/disease_families_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `family` | string | disease family (scripts/disease_families.py), in order of significance. |
| `n_terms` | integer | DisGeNET terms assigned to it. |
| `n_genes` | integer | genes of the universe in at least one of its terms. |
| `n_methylated_genes` | integer | of these, genes of methylated proteins. |

### `supporting/density_windows_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `window_aa` | integer | window width in residues. |
| `start` | integer | first residue of the densest window of that width. |
| `end` | integer | last residue. |
| `sites` | integer | methylarginine sites in it. |
| `arginines` | integer | arginines in it. |
| `null_mean_sites` | number | mean of the most sites in a window under the null (sites placed on the protein's arginines at random). |
| `log2_obs_over_null` | number | log2 of observed over null mean. |
| `p` | number | empirical p value (1,000 permutations). |
| `q` | number | Benjamini-Hochberg q value over all protein-window tests. Only windows with q <= 0.05 are listed. |
| `site_positions` | list&lt;integer&gt; | the sites in the window. |

### `supporting/condensates_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `condensate_id` | string | CD-CODE condensate id. |
| `condensate_name` | string | its name. |
| `condensate_type` | string | its type in CD-CODE. |
| `n_members` | integer | human proteins CD-CODE lists as members. |
| `n_members_methylated` | integer | of these, proteins of the site table. |

### `supporting/hpa_locations_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `gene` | string | gene as HPA gives it. |
| `ensembl_gene` | string | Ensembl gene HPA maps to the protein. |
| `location` | string | subcellular location. |
| `main` | boolean | a main location (otherwise additional). |
| `reliability` | string | HPA reliability: Enhanced, Supported, Approved or Uncertain. |

### `supporting/other_ptm_sites_v1.1.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `position` | integer | position on the reference sequence. |
| `residue` | string | the modified residue. |
| `ptm_type` | code | the modification. Codes: `phospho`, `lys_acyl`, `ub_sumo`, `lys_methyl`. |
| `source` | string | ptm_master: the iPTMnet-format PTM master table (rows on the canonical accession, residue checked against the reference sequence). |
| `evidence` | list&lt;string&gt; | the databases the rows come from, as the master table names them. |
| `n_records` | integer | rows of the master table. |
