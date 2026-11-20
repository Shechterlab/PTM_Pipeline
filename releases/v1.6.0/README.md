# Human methylarginine sites, context marks, site-level interaction evidence, and cancer and germline variants, v1.6.0

v1.6.0 holds everything of v1.5.1 under this release's names, its data unchanged: the methylarginine data of v1.4.3, the tables of context marks and the site-level interaction tables of every site set. It adds, at every candidate residue of the site sets, the somatic changes in TCGA tumors and DepMap cell lines and the germline variants of ClinVar and gnomAD, each with the site sets the residue belongs to, and per-protein background counts. See CHANGES.md for what changed.

## New: somatic changes at every candidate residue

`supporting/tcga_mc3_changes_<residue>_v1.6.0.tsv.gz` and `supporting/depmap_24q4_changes_<residue>_v1.6.0.tsv.gz` list the somatic missense and stop changes at every arginine, lysine, serine, threonine and tyrosine of the reviewed canonical proteins, a known site or not: in TCGA tumors (the MC3 calls of the PanCancer Atlas, 10,294 tumor samples in 33 projects) and in DepMap 24Q4 cell lines (1,929 lines with mutation calls). One row per tumor sample or cell line and change. The column `site_sets` names the site sets of this release with a site at that residue (empty: a candidate residue, not a known site), so a set's changes can be read at every candidate residue or at its known sites only.

| Residue | TCGA changes | at known sites | DepMap changes | at known sites |
|---|---:|---:|---:|---:|
| Arginine (R) | 274,408 | 4,522 | 53,209 | 1,041 |
| Lysine (K) | 61,510 | 7,926 | 20,835 | 2,841 |
| Serine (S) | 142,576 | 14,218 | 44,341 | 4,203 |
| Threonine (T) | 59,012 | 2,810 | 24,605 | 1,262 |
| Tyrosine (Y) | 25,757 | 919 | 13,533 | 493 |

| Site set | Sites | Sites changed in TCGA tumors | Sites changed in DepMap cell lines |
|---|---:|---:|---:|
| Methylarginine (R) | 11,632 | 3,100 | 849 |
| Methyllysine (K) | 5,178 | 415 | 151 |
| Lysine acylation, ubiquitination and SUMOylation (K) | 90,876 | 6,984 | 2,466 |
| Phosphorylation (S, T, Y) | 142,415 | 15,479 | 5,303 |

`supporting/cancer_change_protein_background_v1.6.0.tsv.gz` gives, for each protein and residue class (and over all its residues), how many residues it has and the changes, samples or lines and positions on them, so a site's changes can be compared with those of the other residues of its kind in the same protein. Arginines are changed far more often than the other residues (their codons hold CpG dinucleotides, which mutate often), so such comparisons are best made within a residue class.

The rules are those of the arginine layer of earlier releases (`scripts/44b_extract_tcga_mc3.py`, `scripts/44_stage_depmap.py`): TCGA calls with FILTER PASS (in OV and LAML also `wga` alone), read on a transcript UniProt assigns to the canonical isoform, else on the MAF's own transcript placed through its entry or gene symbol; DepMap calls on the canonical isoform; a change counts only where the 2026_03 canonical sequence has that residue. `scripts/66_stage_cancer_changes.py` stages the K, S, T and Y tables and, rerun on arginines, writes the arginine tables byte for byte. TCGA data are GDC open access (NIH Genomic Data Sharing policy; acknowledge TCGA and the MC3 publication, Ellrott et al., Cell Syst 6:271, 2018, and make no attempt to identify participants); DepMap 24Q4 is CC BY 4.0 (doi:10.25452/figshare.plus.27993248).

`supporting/tcga_mc3_projects_v1.6.0.tsv` names the TCGA projects and their tumor samples, and `supporting/depmap_24q4_change_models_v1.6.0.tsv` the cell lines (name, lineage, primary disease).

## New: germline variants at every candidate residue (ClinVar, gnomAD)

`supporting/clinvar_changes_<residue>_v1.6.0.tsv.gz` and `supporting/gnomad_v4_1_changes_<residue>_v1.6.0.tsv.gz` list the human germline variants that change an arginine, lysine, serine, threonine or tyrosine of a reviewed canonical protein (a substitution or a stop), a known site or not, with `site_sets` as in the cancer tables. ClinVar (files of 2026-11-15 (hgvs4variation), 2026-11-16 (variant_summary)): one row per variation, with its germline classification, review status and stars, submitters, origin, conditions, somatic clinical impact and oncogenicity; the variant counts where one of its protein changes is written on a RefSeq protein UniProt assigns to the entry's displayed isoform. gnomAD v4.1 (730,947 exomes): one row per variant that passes every gnomAD filter, with its allele count and number, frequency, homozygotes, the ancestry group with the highest frequency and the highest filtering allele frequency; the change is VEP's, on a transcript UniProt assigns to the displayed isoform. Either way a variant counts only where the 2026_03 canonical sequence has that residue.

| Residue | ClinVar variants | at known sites | of which pathogenic or likely pathogenic | gnomAD variants | at known sites |
|---|---:|---:|---:|---:|---:|
| Arginine (R) | 335,933 | 6,026 | 412 | 1,133,321 | 20,550 |
| Lysine (K) | 99,522 | 15,986 | 995 | 517,573 | 75,792 |
| Serine (S) | 186,004 | 23,071 | 947 | 965,320 | 108,765 |
| Threonine (T) | 142,865 | 9,048 | 194 | 678,065 | 37,713 |
| Tyrosine (Y) | 64,910 | 2,645 | 522 | 307,801 | 10,081 |

| Site set | Sites | Sites with a ClinVar variant | with a pathogenic or likely pathogenic one | Sites with a gnomAD variant |
|---|---:|---:|---:|---:|
| Methylarginine (R) | 11,632 | 4,076 | 328 | 9,139 |
| Methyllysine (K) | 5,178 | 693 | 46 | 2,819 |
| Lysine acylation, ubiquitination and SUMOylation (K) | 90,876 | 12,028 | 866 | 49,864 |
| Phosphorylation (S, T, Y) | 142,415 | 27,085 | 1,316 | 96,387 |

`supporting/clinvar_protein_background_v1.6.0.tsv.gz` and `supporting/gnomad_protein_background_v1.6.0.tsv.gz` give, for each protein and residue class, the residues and the variants on them (and, for ClinVar, those pathogenic or likely pathogenic; for gnomAD, the summed allele count).

NCBI places no restrictions on the use or distribution of ClinVar data (submitters may claim rights in what they submitted); cite ClinVar (Landrum et al., Nucleic Acids Res 46:D1062, 2018). The gnomAD primary data are CC0; cite gnomAD v4 (Chen et al., Nature 625:92, 2024). gnomAD's SpliceAI and CADD annotations, which carry non-commercial terms, are not used. gnomAD was read from Google's public copy (`gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/`), streamed (`scripts/69_stage_gnomad.py`); ClinVar from its FTP tables (`scripts/67_stage_clinvar.py`).

## Counts

Sites are the distinct residues of each site set's tables. High-confidence sites are the website's control, not a tier of the release: the release tables keep every row. For the context marks a site is high confidence when at least two independent studies report it from experimental records: each PubMed ID counts once, a direct dataset (Berryhill 2023, Ochoa 2020, Hendriks 2017, Akimov 2018) counts as one study, and so as the paper a database cites, and UniProt and iPTMnet count only through the papers they cite; a residue's studies are pooled over all its rows of one modification. For the lysine set the studies must report the same modification (each acyl type, ubiquitin and SUMO apart), so a lysine that one study reports ubiquitinated and another SUMOylated is not high confidence; methylation degrees (Kme1, Kme2, Kme3) are one modification. For methylarginine it is the existing control: two or more studies, a laboratory's studies counted once, PhosphoSitePlus, UniProt, IEDB and ProMetheusDB not counted. Sites with interaction evidence have a row of their set's `site_interactions` table that counts as evidence; the methylarginine tables hold more kinds of evidence (predicted interfaces and contacts, IntAct, RBDmap, pCLAP), so their count from only the kinds the other sets have is given too.

| Site set | Sites | High-confidence sites | Sites with interaction evidence |
|---|---:|---:|---:|
| Methylarginine (R) | 11,632 | 2,985 | 4,530 (3,069 from the kinds of evidence of the other sets) |
| Methyllysine (K) | 5,178 | 83 | 485 |
| Lysine acylation, ubiquitination and SUMOylation (K) | 90,876 | 1,938 | 10,878 |
| &nbsp;&nbsp;of which acetylation | 4,597 | 215 | 1,061 |
| &nbsp;&nbsp;of which lysine acylation (all acyl types) | 5,504 | 218 | 1,294 |
| &nbsp;&nbsp;of which ubiquitination | 63,480 | 375 | 8,150 |
| &nbsp;&nbsp;of which SUMOylation | 40,249 | 1,370 | 5,201 |
| Phosphorylation (S, T, Y) | 142,415 | 42,804 | 10,673 |

## Site-level interaction evidence

The tables are made by `scripts/64_export_site_interactions.py`, which runs the interface, crosslink and structure-contact steps of the methylarginine pipeline (`scripts/37_export_interaction_data.py`) with the modified residue as a parameter (`config/site_sets.json`). Run on the methylarginines it gives the methylarginine release tables (`tests/test_site_interactions.py`).

### Site sets

| Site set | Residues | Sites from | Sites | with a PDB interface | with RNA or DNA within 4 Å in a structure | with evidence from cells |
|---|---|---|---|---|---|---|
| methylarginine | R | `human_methylarginine_sites_v1.6.0.tsv` | 11,632 | 1,396 | 344 | 4,330 |
| methyllysine | K | `human_methyllysine_sites_v1.6.0.tsv` | 5,178 | 419 | 86 | 114 |
| lysine | K | `human_acetylation_sites_v1.6.0.tsv`; `human_lysine_acylation_sites_v1.6.0.tsv`; `human_ubiquitination_sites_v1.6.0.tsv`; `human_sumoylation_sites_v1.6.0.tsv` | 90,876 | 9,529 | 1,706 | 2,429 |
| phospho | S, T, Y | `human_phosphorylation_sites_v1.6.0.tsv` | 142,415 | 6,635 | 585 | 5,042 |

### Tables

One file per site set and table, `supporting/<site set>_<table>_v1.6.0.tsv`, in the columns of the methylarginine tables of the same name (`releases/v1.5.0/data_dictionary.tsv`). A site is `<accession>:<residue><position>`.

| Table | One row per |
|---|---|
| `site_interactions` | piece of evidence at a site: a PDB interface with a protein, RNA, DNA or another copy of the protein (PDBe-KB), an RNA or DNA crosslink, a CAPRI peptide, an RBR-ID residue, a UniProt DNA-binding region |
| `site_interaction_records` | site and record: each PDB entry or study behind `site_interactions`, with its method, paper and what the entry models at the site |
| `nucleic_acid_contacts` | site and RNA or DNA chain within 4 Å of its side chain in a PDB entry: the contacted nucleotides and the bound sequence around them |
| `nucleic_acid_proximal_sites` | site and complex with RNA or DNA within 8 Å of the side chain and none within 4 Å |
| `nucleic_acid_near_sites` | site and complex in which a residue at i ± 1–5, not the site, is within 4 Å of RNA or DNA (the side chain not within 8 Å) |
| `residue_evidence` | run of residues of a protein of the set with the same evidence |

Side-chain distances are measured from the atoms that carry the modification and the carbons next to them: lysine CD, CE, NZ; serine CB, OG; threonine CB, OG1, CG2; tyrosine CE1, CE2, CZ, OH.

#### record_class

| Code | Meaning |
|---|---|
| `binds_modified` | binds with the site modified in the experiment: a PDB structure modeling a modified residue there (the components of config/site_sets.json modified_components) |
| `binds_unmodified` | binds with the unmodified residue at the site: a PDB structure modeling it there |
| `binds_other_residue` | binds with another residue at the site: a PDB structure of a mutated or variant construct (state_detail names it) |
| `binds_state_not_read` | a PDB interface at the site whose residue there could not be read (no entity of the entry aligns to the site) |
| `nucleic_acid` | RNA or DNA: contact in a PDB structure, a crosslink or RNA-bound peptide in cells, an RBR-ID residue, or a UniProt DNA-binding region (annotation) |
| `self` | another copy of the same protein in a PDB structure |

#### site_state

| Code | Meaning |
|---|---|
| `lysine` | an unmodified lysine (PDB residue LYS) |
| `serine` | an unmodified serine (PDB residue SER) |
| `threonine` | an unmodified threonine (PDB residue THR) |
| `tyrosine` | an unmodified tyrosine (PDB residue TYR) |
| `Kac` | lysine modified, state Kac (PDB component ALY) |
| `Kme1` | lysine modified, state Kme1 (PDB component MLZ) |
| `Kme2` | lysine modified, state Kme2 (PDB component MLY) |
| `Kme3` | lysine modified, state Kme3 (PDB component M3L) |
| `Sph` | serine modified, state Sph (PDB component SEP) |
| `Tph` | threonine modified, state Tph (PDB component TPO) |
| `Yph` | tyrosine modified, state Yph (PDB component PTR) |
| `other_residue` | another residue in the construct (state_detail names it) |
| `not_read` | a structure whose residue at the site could not be read |
| `not_recorded` | not recorded: a crosslink, RNA-bound peptide or annotation |

The tables of methylarginine hold more kinds of evidence, classes and states than these (predicted interfaces and contacts, IntAct features, RBDmap and pCLAP): `data_dictionary.tsv` and the methylarginine section below describe them. In the table above, their evidence from cells counts those kinds too.

### Not in the tables of methyllysine, lysine, phospho

- RBDmap fragments and pCLAP peptides: `scripts/64` does not read them for these sets. The methylarginine tables carry them, with their extents and sequences: RBDmap from Castello et al., Mol Cell 63:696 (2016), doi:10.1016/j.molcel.2016.06.029 (CC BY 4.0); pCLAP, the RNA-binding peptides of Mullari et al., J Proteome Res 16:2762 (2017), doi:10.1021/acs.jproteome.7b00042, called by the pipeline from PRIDE PXD005739 with the rule of Mullari et al., Nat Commun 14:4348 (2023).
- Predicted interfaces and contacts (ProtVar, AlphaFold DB) and IntAct features: the pipeline computes them for methylarginines only.
- Enzyme-substrate records (kinases, writers and erasers).
- Statistics: no test is run on these tables.

### Sources and licenses

| Source | Gives | License |
|---|---|---|
| PDBe-KB interface residues (PDBe-KB consortium, Nucleic Acids Res 50:D534, 2022) | the residues at protein, RNA and DNA interfaces in PDB entries | EMBL-EBI terms of use (free use with attribution) |
| SIFTS (Dana et al., Nucleic Acids Res 47:D482, 2019) | UniProt numbering of PDB chains | EMBL-EBI terms of use |
| PDB, read through the RCSB ModelServer and Data API and the PDBe API (Burley et al., Nucleic Acids Res 51:D488, 2023) | side-chain contacts, bound sequences, entry titles, methods and citations, the residue each entry models at a site | CC0 1.0 |
| pepR-MS (Sha, Kuster & Trendel, bioRxiv 2025, doi:10.1101/2025.09.05.674461), a preprint that re-analyses Kramer 2014, XRNAX (Trendel 2019), RBS-ID (Bae 2020) and pRBS-ID (Bae 2021) | RNA-crosslinked residues | CC BY-NC |
| CAPRI (Panhale et al., Nat Commun 10:2682, 2019) | RNA-crosslinked peptides and the peptides next to them | CC BY 4.0 |
| RBR-ID (Oksuz et al., Mol Cell 83:2449, 2023; PRIDE) | residues under RBR-ID+ peptides, regions of high residue score | CC0 |
| Trendel et al., Cell 188:4424 (2025) | DNA-crosslinked residues and peptides | CC BY-NC 4.0 |
| UniProt (UniProt Consortium, Nucleic Acids Res 53:D609, 2025) | DNA-binding regions; sequences and gene names | CC BY 4.0 |

The tables are under CC BY-NC-SA 4.0, as the release. Cite the sources of the evidence used.

### Files

| File | Rows |
|---|---|
| `supporting/lysine_nucleic_acid_contacts_v1.6.0.tsv` | 28,471 |
| `supporting/lysine_nucleic_acid_near_sites_v1.6.0.tsv` | 7,652 |
| `supporting/lysine_nucleic_acid_proximal_sites_v1.6.0.tsv` | 10,820 |
| `supporting/lysine_residue_evidence_v1.6.0.tsv` | 179,341 |
| `supporting/lysine_site_interaction_records_v1.6.0.tsv` | 118,646 |
| `supporting/lysine_site_interactions_v1.6.0.tsv` | 19,011 |
| `supporting/methylarginine_nucleic_acid_contacts_v1.6.0.tsv` | 14,212 |
| `supporting/methylarginine_nucleic_acid_near_sites_v1.6.0.tsv` | 943 |
| `supporting/methylarginine_nucleic_acid_proximal_sites_v1.6.0.tsv` | 1,519 |
| `supporting/methylarginine_residue_evidence_v1.6.0.tsv` | 222,880 |
| `supporting/methylarginine_site_interaction_records_v1.6.0.tsv` | 52,060 |
| `supporting/methylarginine_site_interactions_v1.6.0.tsv` | 22,899 |
| `supporting/methyllysine_nucleic_acid_contacts_v1.6.0.tsv` | 1,466 |
| `supporting/methyllysine_nucleic_acid_near_sites_v1.6.0.tsv` | 678 |
| `supporting/methyllysine_nucleic_acid_proximal_sites_v1.6.0.tsv` | 1,326 |
| `supporting/methyllysine_residue_evidence_v1.6.0.tsv` | 58,802 |
| `supporting/methyllysine_site_interaction_records_v1.6.0.tsv` | 8,367 |
| `supporting/methyllysine_site_interactions_v1.6.0.tsv` | 1,358 |
| `supporting/phospho_nucleic_acid_contacts_v1.6.0.tsv` | 14,030 |
| `supporting/phospho_nucleic_acid_near_sites_v1.6.0.tsv` | 4,093 |
| `supporting/phospho_nucleic_acid_proximal_sites_v1.6.0.tsv` | 4,954 |
| `supporting/phospho_residue_evidence_v1.6.0.tsv` | 198,172 |
| `supporting/phospho_site_interaction_records_v1.6.0.tsv` | 81,420 |
| `supporting/phospho_site_interactions_v1.6.0.tsv` | 18,389 |

## Context marks

| File | One row per |
|---|---|
| `human_methyllysine_sites_v1.6.0.tsv` | methyllysine residue and state (5,746 rows on 5,178 residues of 2,843 proteins) |
| `human_phosphorylation_sites_v1.6.0.tsv` | phosphorylation residue and state (142,415 rows on 142,415 residues of 12,588 proteins) |
| `human_lysine_acylation_sites_v1.6.0.tsv` | lysine acylation residue and acyl type (8,008 rows on 5,504 residues of 1,945 proteins) |
| `human_acetylation_sites_v1.6.0.tsv` | acetylation residue and acyl type (4,597 rows on 4,597 residues of 1,834 proteins) |
| `human_ubiquitination_sites_v1.6.0.tsv` | ubiquitination residue and state (63,480 rows on 63,480 residues of 9,032 proteins) |
| `human_sumoylation_sites_v1.6.0.tsv` | sumoylation residue and state (40,249 rows on 40,249 residues of 6,410 proteins) |
| `human_methylarginine_crosscheck_v1.6.0.tsv` | arginine methylation record of UniProt or iPTMnet (2,471 rows), with whether the methylarginine data has the arginine. A cross-check, not a second methylarginine table |

Columns are in `data_dictionary.tsv`. Each row says which sources report it (`sources`), the best `evidence_class`, UniProt's ECO codes, iPTMnet's upstream databases (`iptmnet_sources`), the PubMed IDs and, in its own columns, the figures of Berryhill 2023 (`berryhill_*`), of Ochoa 2020 (`functional_score`, `best_localization_probability`, `best_pep`, `biological_samples`, `spectral_counts`), of Hendriks 2017 (`hendriks_*`, sumoylation) and of Akimov 2018 (`akimov_*`, ubiquitination) and, in the methyllysine, phosphorylation, acylation, ubiquitination and sumoylation tables, iPTMnet's own score (`iptmnet_score`); each of these is empty on the rows of other sources. A site can hide `by-similarity` (UniProt ECO:0000250 and other non-experimental records) and `unspecified` rows by filtering on `evidence_class`.

### Sources and licenses

| Source | Release | License | Terms |
|---|---|---|---|
| UniProt | 2026_04 | CC BY 4.0 | https://www.uniprot.org/help/license |
| iPTMnet | 6.2 | CC BY-NC-SA 4.0 | https://research.bioinformatics.udel.edu/iptmnet/license |
| Berryhill 2023 | Sci Rep 13:377 (2023) | article CC BY 4.0; supplements carry no separate notice | https://www.nature.com/articles/s41598-022-27175-x |
| Ochoa 2020 | Nat Biotechnol 38:365 (2020) | article not open access (exclusive license to Springer Nature America); supplements carry no license | https://www.nature.com/articles/s41587-019-0344-3 |
| Hendriks 2017 | Nat Struct Mol Biol 24:325 (2017) | article not open access (copyright 2017 Springer Nature America); supplements carry no notice | https://www.nature.com/articles/nsmb.3366 |
| Akimov 2018 | Nat Struct Mol Biol 25:631 (2018) | article not open access (copyright 2018 The Author(s)); supplements carry no notice | https://www.nature.com/articles/s41594-018-0084-y |

The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)). UniProt (CC BY 4.0) is compatible with it. iPTMnet's license page says CC BY-NC-SA 4.0; its download readme also mentions CC BY 4.0 in its last line. The release follows the stricter terms: share-alike and non-commercial. iPTMnet asks that the upstream databases of its `source` column are cited when its data are redistributed: each row keeps them (`iptmnet_sources`), and they are HPRD, UniProt, RLIMS-P text mining, SIGNOR, neXtProt, IEDB, Protein Ontology and IntAct (terms are theirs). Not used: PhosphoSitePlus, dbPTM, qPTM and CPLM 4.0. CPLM 4.0 states no license (cplm.biocuckoo.cn carries a copyright notice, "All Rights Reserved") and compiles sites from other databases as well as from papers, so it is left out.

**Berryhill 2023** (Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. "Global lysine methylome profiling using systematically characterized affinity reagents", Sci Rep 13:377, 2023, doi:10.1038/s41598-022-27175-x). The article is CC BY 4.0; its supplements carry no separate license notice. Its Supplementary Table 5 is included, cited. It contributes 5,072 methyllysine rows (states me1, me2, me3 kept; `evidence_class` experimental, as the sites are identified by mass spectrometry in HEK293T or U2OS cells; 4,018 with a site probability) on 2,743 proteins. Of the 6,576 rows of the table, 84 sites were moved to the lysine of their peptide on the canonical sequence and 17 rows were rejected (8 accession not in reference, 3 peptide not in canonical, 6 peptide not unique). The other sheets of the supplements, which carry third-party information, are not used.

**Ochoa 2020** (Ochoa D, Jarnuczak AF, Viéitez C, et al. "The functional landscape of the human phosphoproteome", Nat Biotechnol 38:365-373, 2020, doi:10.1038/s41587-019-0344-3). The article is **not open access** ("© 2019 The Author(s), under exclusive licence to Springer Nature America, Inc."), and its supplements carry no license of their own. The derived table (accession, position, residue, functional score, best localization probability, best PEP, biological samples, spectral counts) is published here, cited. The CC BY-NC-SA 4.0 of this release does not grant rights in Ochoa et al.'s data beyond what Springer Nature and the authors allow; a reader who needs the phosphorylation table without it removes that source from `sources`, and the rows that only it reports. It contributes 112,405 phosphorylation rows (state `phospho`; `evidence_class` experimental, mass spectrometry) on 10,275 proteins; of its 116,258 sites 3,853 were rejected on the reference check (1,200 accession not in reference, 105 bad position, 2,548 residue mismatch; no window is given, so none was relocated).

**Hendriks 2017** (Hendriks IA, Lyon D, Young C, Jensen LJ, Vertegaal ACO, Nielsen ML. "Site-specific mapping of the human SUMO proteome reveals co-modification with phosphorylation", Nat Struct Mol Biol 24:325-336, 2017, doi:10.1038/nsmb.3366; PRIDE PXD005296). The article is not open access (copyright "2017 Springer Nature America, Inc."), and its supplementary table carries no notice of its own. Its Supplementary Table 1 (all SUMO-2 acceptor lysines identified in HeLa and U2OS cells, as filtered by the paper) is included, cited, with the localization probability (`hendriks_localization_probability`) and the cell lines with an MS/MS detection (`hendriks_cell_lines`). Of the 40,765 sites of the table, 39,935 were placed on the canonical sequence of the reference (38,955 kept, 529 relocated, 451 secondary accession; kept: a protein of the site's protein group is a canonical entry with the lysine at the given position and the +/-7 residue window agreeing; secondary accession: the same after the current accession replaces a merged one; relocated: the window occurs once in the canonical sequence of a protein of the group, mostly sites given on an isoform) and 830 were left out (329 no canonical accession, 496 window not in canonical, 5 window not unique); 77 placed sites fell on a lysine another of its sites also gives. It contributes 39,858 sumoylation rows (state `sumo`; `evidence_class` experimental, mass spectrometry) on 6,313 proteins: 5,149 residues that other sources also report and 34,709 that only it reports.

**Akimov 2018** (Akimov V, Barrio-Hernandez I, Hansen SVF, et al. "UbiSite approach for comprehensive mapping of lysine and N-terminal ubiquitination sites", Nat Struct Mol Biol 25:631-640, 2018, doi:10.1038/s41594-018-0084-y). The article is not open access (copyright "2018 The Author(s)"), and its supplementary table carries no notice of its own. Its Supplementary Table 1 (the lysine ubiquitination sites of the UbiSite-enriched ubiquitinomes of Hep2 and Jurkat cells, both digests combined, as filtered by the paper) is included, cited, with the best localization probability (`akimov_localization_probability`) and the cell lines in which the site has an intensity value (`akimov_cell_lines`). The N-terminal ubiquitination sites of the paper (its Supplementary Table 6) are not included: the table holds lysine sites only. Of the 63,489 sites of the table, 62,951 were placed on the canonical sequence of the reference (61,611 kept, 1,332 relocated, 8 secondary accession; kept: a protein of the site's protein group is a canonical entry with the lysine at the given position and the +/-7 residue window agreeing; secondary accession: the same after the current accession replaces a merged one; relocated: the window occurs once in the canonical sequence of a protein of the group, mostly sites given on an isoform) and 538 were left out (3 bad window, 36 no canonical accession, 495 window not in canonical, 4 window not unique); 80 placed sites fell on a lysine another of its sites also gives. It contributes 62,871 ubiquitination rows (state `ub`; `evidence_class` experimental, mass spectrometry) on 8,927 proteins: 349 residues that other sources also report and 62,522 that only it reports.


**iPTMnet score** (`iptmnet_score`). iPTMnet's download has a `score.txt` (iPTMnet 6.2, https://research.bioinformatics.udel.edu/iptmnet_data/files/current/score.txt); its readme says only "The score is calculated based on all the evidence sources as shown on the website". It is an integer from 0 to 4 and is iPTMnet's own: this release carries it unchanged and does not explain it. `score.txt` scores the site (the rows without an enzyme) under the plain accession and, less often, the accession with the suffix -1 (its isoform 1); the plain accession's score is used. A row has the score when its `sources` include iPTMnet and `score.txt` has a row for the residue: methyllysine 231 of 234; phosphorylation 63,586 of 64,117; lysine acylation 2,657 of 2,672; ubiquitination 609 of 610; sumoylation 33 of 33.

### Counts

Rows by state, evidence class and the sources that report them (iPTMnet 6.2 gives no methylation degree, so its rows are state `me` and do not meet UniProt's `me1`, `me2` and `me3` rows; ubiquitination rows are state `ub` and sumoylation rows state `sumo`, UniProt's SUMO paralog not kept).

**methyllysine** (5,746 rows). Residues by source, any state: Berryhill 2023 only 4,766, UniProt + Berryhill 2023 15, UniProt + iPTMnet 175, UniProt + iPTMnet + Berryhill 2023 44, UniProt only 163, iPTMnet only 15.

| State | experimental | by-similarity | unspecified | UniProt only | iPTMnet only | Berryhill 2023 only | UniProt + iPTMnet | UniProt + Berryhill 2023 | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| me | 27 | 18 | 206 | 17 | 223 | 0 | 11 | 0 | 251 |
| me1 | 2,230 | 72 | 0 | 218 | 0 | 2,046 | 0 | 38 | 2,302 |
| me2 | 1,348 | 36 | 0 | 102 | 0 | 1,257 | 0 | 25 | 1,384 |
| me3 | 1,754 | 55 | 0 | 103 | 0 | 1,687 | 0 | 19 | 1,809 |
| all | 5,359 | 181 | 206 | 440 | 223 | 4,990 | 11 | 82 | 5,746 |

**phosphorylation** (142,415 rows). Residues by source, any state: Ochoa 2020 only 70,281, UniProt + Ochoa 2020 3,461, UniProt + iPTMnet 7,279, UniProt + iPTMnet + Ochoa 2020 25,632, UniProt only 4,556, iPTMnet + Ochoa 2020 13,031, iPTMnet only 18,175.

| State | experimental | text-mined | by-similarity | unspecified | UniProt only | iPTMnet only | Ochoa 2020 only | UniProt + iPTMnet | UniProt + Ochoa 2020 | iPTMnet + Ochoa 2020 | UniProt + iPTMnet + Ochoa 2020 | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| phospho | 135,936 | 2,018 | 4,393 | 68 | 4,556 | 18,175 | 70,281 | 7,279 | 3,461 | 13,031 | 25,632 | 142,415 |
| all | 135,936 | 2,018 | 4,393 | 68 | 4,556 | 18,175 | 70,281 | 7,279 | 3,461 | 13,031 | 25,632 | 142,415 |

**lysine_acylation** (8,008 rows). Residues by source, any state: UniProt + iPTMnet 2,475, UniProt only 2,832, iPTMnet only 197.

| State | experimental | by-similarity | unspecified | UniProt only | iPTMnet only | UniProt + iPTMnet | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2-hydroxyisobutyryl | 365 | 54 | 0 | 419 | 0 | 0 | 419 |
| acetyl | 2,710 | 1,882 | 5 | 1,925 | 242 | 2,430 | 4,597 |
| benzoyl | 115 | 45 | 0 | 160 | 0 | 0 | 160 |
| beta-hydroxybutyryl | 221 | 57 | 0 | 278 | 0 | 0 | 278 |
| butyryl | 50 | 12 | 0 | 62 | 0 | 0 | 62 |
| crotonyl | 210 | 14 | 0 | 224 | 0 | 0 | 224 |
| glutaryl | 227 | 37 | 0 | 264 | 0 | 0 | 264 |
| isonicotinyl | 145 | 44 | 0 | 189 | 0 | 0 | 189 |
| lactoyl | 197 | 75 | 0 | 272 | 0 | 0 | 272 |
| malonyl | 39 | 6 | 0 | 45 | 0 | 0 | 45 |
| methacryl | 160 | 55 | 0 | 215 | 0 | 0 | 215 |
| propionyl | 10 | 9 | 0 | 19 | 0 | 0 | 19 |
| succinyl | 84 | 1,180 | 0 | 1,264 | 0 | 0 | 1,264 |
| all | 4,533 | 3,470 | 5 | 5,336 | 242 | 2,430 | 8,008 |

**ubiquitination** (63,480 rows). Residues by source, any state: Akimov 2018 only 62,522, UniProt + Akimov 2018 103, UniProt + iPTMnet 358, UniProt + iPTMnet + Akimov 2018 245, UniProt only 245, iPTMnet + Akimov 2018 1, iPTMnet only 6.

| State | experimental | by-similarity | UniProt only | iPTMnet only | Akimov 2018 only | UniProt + iPTMnet | UniProt + Akimov 2018 | iPTMnet + Akimov 2018 | UniProt + iPTMnet + Akimov 2018 | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ub | 63,299 | 181 | 245 | 6 | 62,522 | 358 | 103 | 1 | 245 | 63,480 |
| all | 63,299 | 181 | 245 | 6 | 62,522 | 358 | 103 | 1 | 245 | 63,480 |

**sumoylation** (40,249 rows). Residues by source, any state: Hendriks 2017 only 34,709, UniProt + Hendriks 2017 5,128, UniProt + iPTMnet 6, UniProt + iPTMnet + Hendriks 2017 18, UniProt only 379, iPTMnet + Hendriks 2017 3, iPTMnet only 6.

| State | experimental | by-similarity | UniProt only | iPTMnet only | Hendriks 2017 only | UniProt + iPTMnet | UniProt + Hendriks 2017 | iPTMnet + Hendriks 2017 | UniProt + iPTMnet + Hendriks 2017 | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sumo | 40,090 | 159 | 379 | 6 | 34,709 | 6 | 5,128 | 3 | 18 | 40,249 |
| all | 40,090 | 159 | 379 | 6 | 34,709 | 6 | 5,128 | 3 | 18 | 40,249 |

**methylarginine cross-check** (2,471 rows). Residues by source, any state: UniProt + iPTMnet 742, UniProt only 740, iPTMnet only 65.

| State | experimental | by-similarity | unspecified | UniProt only | iPTMnet only | UniProt + iPTMnet | Rows |
|---|---:|---:|---:|---:|---:|---:|---:|
| me | 87 | 1 | 719 | 0 | 802 | 5 | 807 |
| me1 | 502 | 483 | 0 | 985 | 0 | 0 | 985 |
| me2 | 38 | 33 | 0 | 71 | 0 | 0 | 71 |
| me2a | 308 | 245 | 0 | 553 | 0 | 0 | 553 |
| me2s | 21 | 34 | 0 | 55 | 0 | 0 | 55 |
| all | 956 | 796 | 719 | 1,664 | 802 | 5 | 2,471 |

Arginines of the cross-check in the methylarginine data: lower_confidence 63, no 308, site 1,176.

iPTMnet carries a copy of UniProt: 217 methyllysine, 11,692 phosphorylation, 2,257 lysine acylation, 549 ubiquitination, 0 sumoylation rows have `unip` as their only upstream database, so "UniProt + iPTMnet" counts them as agreement of two sources when they are one.

### Dropped in staging

Each count is of source records. UniProt: 763 on an isoform, 2 residue mismatch; the entries' sequences were compared with the reference (UniProt 2026_04 against 2026_03). iPTMnet: 273,145 not human, 283 not in reference, 18,652 on an isoform, 32,107 other ptm type, 8,890 other residue, 483 residue mismatch. Berryhill 2023: 6,576 rows of Supplementary Table 5, 8 accession not in reference, 3 peptide not in canonical, 6 peptide not unique (and 84 sites relocated by their peptide). Ochoa 2020: 116,258 sites, 1,200 accession not in reference, 105 bad position, 2,548 residue mismatch (nothing relocated). Hendriks 2017: 40,765 sites of Supplementary Table 1, 329 no canonical accession, 496 window not in canonical, 5 window not unique (and 529 relocated by their window). Akimov 2018: 63,489 sites of Supplementary Table 1, 3 bad window, 36 no canonical accession, 495 window not in canonical, 4 window not unique (and 1,332 relocated by their window). iPTMnet scores: 851,361 rows of score.txt, 23,643 enzyme row, 331,220 not in reference, 27,456 on an isoform, 70,695 other ptm type or residue, 2,177 residue mismatch, 517 unparsed site.

### How the context marks were made

`python scripts/54_stage_uniprot_ptm_sites.py`, `python scripts/55_stage_iptmnet_ptm_sites.py`, `python scripts/57_stage_berryhill_lys_methylation.py`, `python scripts/58_stage_ochoa_phosphosites.py`, `python scripts/59_stage_iptmnet_scores.py`, `python scripts/60_stage_hendriks2017_sumo_sites.py` and `python scripts/61_stage_akimov2018_ubiquitination_sites.py` (inputs in `inputs/marks/`, downloads listed in `inputs/marks/raw_files.tsv`), then `python scripts/56_build_context_marks.py --version v1.5.0 --release-date 2026-11-17` at commit `125f10960e39445db5a32ec4a30f5bb8cfed9f29` of https://github.com/Shechterlab/PTM_Pipeline.

### Cite

The paper below, and the sources you use:

- UniProt: UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53:D609-D617 (2025). doi:10.1093/nar/gkae1010
- iPTMnet: Huang H, Arighi CN, Ross KE, et al. iPTMnet: an integrated resource for protein post-translational modification network discovery. Nucleic Acids Res 46:D542-D550 (2018). doi:10.1093/nar/gkx1104
- Berryhill 2023: Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. Global lysine methylome profiling using systematically characterized affinity reagents. Sci Rep 13:377 (2023). doi:10.1038/s41598-022-27175-x
- Ochoa 2020: Ochoa D, Jarnuczak AF, Viéitez C, et al. The functional landscape of the human phosphoproteome. Nat Biotechnol 38:365-373 (2020). doi:10.1038/s41587-019-0344-3
- Hendriks 2017: Hendriks IA, Lyon D, Young C, Jensen LJ, Vertegaal ACO, Nielsen ML. Site-specific mapping of the human SUMO proteome reveals co-modification with phosphorylation. Nat Struct Mol Biol 24:325-336 (2017). doi:10.1038/nsmb.3366
- Akimov 2018: Akimov V, Barrio-Hernandez I, Hansen SVF, et al. UbiSite approach for comprehensive mapping of lysine and N-terminal ubiquitination sites. Nat Struct Mol Biol 25:631-640 (2018). doi:10.1038/s41594-018-0084-y

## Methylarginine data (unchanged from v1.4.3)

What follows is the README of v1.4.3, with its headings one level down. The files it names are in this folder with `_v1.6.0` in place of `_v1.4.3`.

## Human methylarginine sites, v1.4.3

11,632 sites on 4,439 proteins, 11,632 rows: the human protein arginine methylation sites of DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with corrections to how some source records were placed on UniProt sequences (v1.1.0); from v1.2.0 every study is read from its own tables and a site is kept only where its methyl group is localized to one arginine and the sites of 6 more studies (Hartel 2019, Geoghegan 2015, Sylvestersen 2014, Lu 2025, Ma 2023, Lu 2024). [CHANGES.md](CHANGES.md) says what changed since v1.4.2 and why.

Every site has its sequence, methylation evidence, structural context, the studies that report it, and the source records behind it. Positions are on UniProt release 2026_03 (reviewed human reference proteome UP000005640).

### Files

| File | One row per |
|---|---|
| `human_methylarginine_sites_v1.6.0.tsv` | site, including the 0 sites kept pending (flag `one_sided_window`) |
| `human_methylarginine_sites_v1.6.0.xlsx` | the same table as a spreadsheet |
| `human_methylarginine_proteins_v1.6.0.tsv` | protein: names, site counts, functional classes, arginine enrichment, site density and clusters, condensates, locations, disease families, other PTMs, interactions |
| `human_methylarginine_lower_confidence_sites_v1.6.0.tsv` | lower-confidence site: an arginine a study reports as methylated, or an earlier release had as a site, whose records each fail a site test, with the results of the tests |
| `studies_v1.6.0.tsv` | study or database: citation, DOI, terms, counts |
| `statistics_v1.6.0.tsv` | enrichment test behind the summaries: classes, motifs, disorder, domains, disease families, condensates, locations, interactions |
| `statistics_with_lower_confidence_v1.6.0.tsv` | the same tests, run again with the lower-confidence sites whose own evidence covers the arginine and leaves the methyl group on an arginine |
| `supporting/source_records_v1.6.0.tsv` | source record: what the study reported, its evidence, where each version placed it |
| `supporting/lower_confidence_evidence_v1.6.0.tsv` | source record that adds no site (not localized, or not on one site of a protein it names), with where its methyl group could be |
| `supporting/site_changes_v1.5.0_to_v1.6.0.tsv` | site of either version: kept, removed or added, and why |
| `supporting/published_sites_v1.6.0.tsv` | site of v1.0.0 or v1.1.0: what it is in this release (a site, a lower-confidence site, or not on the reference) |
| `supporting/maron2021_detection_v1.6.0.tsv` | site detected in Maron et al. 2021 Table S2 or S3 (per antibody) |
| `supporting/maron2021_inhibitors_v1.6.0.tsv` | site, antibody and PRMT inhibitor in Maron et al. 2021 Table S3 |
| `supporting/carm1_ko_v1.6.0.tsv` | site and cell line in the CARM1 knockout data of Shishkova et al. 2017 |
| `supporting/methylarginine_site_interactions_v1.6.0.tsv` | piece of interaction evidence at a site: crosslinks and RNA-bound peptides (RBDmap, CAPRI), PDB interfaces, predicted interfaces and contacts, IntAct records |
| `supporting/methylarginine_residue_evidence_v1.6.0.tsv` | run of residues of a protein with the same interaction evidence |
| `supporting/protein_partners_v1.6.0.tsv` | protein and STRING physical partner (score >= 0.7) |
| `supporting/complexes_v1.6.0.tsv` | Complex Portal or hu.MAP complex with a methylated member |
| `supporting/methylarginine_nucleic_acid_contacts_v1.6.0.tsv` | site and RNA or DNA chain within 4 Å of its side chain in a PDB structure: the contacted nucleotides and the bound sequence around them |
| `supporting/arginine_nucleic_acid_distances_v1.6.0.tsv` | arginine of a methylated protein modeled in a PDB entry with RNA or DNA: its minimum side-chain distance to each, methylated or not |
| `supporting/methylarginine_nucleic_acid_proximal_sites_v1.6.0.tsv` | site, complex and kind of nucleic acid in which the arginine side chain is within 8 Å, not 4 Å, of RNA or DNA: the bound sequence around the closest nucleotide |
| `supporting/methylarginine_nucleic_acid_near_sites_v1.6.0.tsv` | site and complex in which a residue at i ± 1–5, not the arginine side chain, is within 4 Å of RNA or DNA (the side chain not within 8 Å): the bound sequence there |
| `supporting/nucleic_acid_regions_v1.6.0.tsv` | run of residues of a protein within 4 Å of RNA or DNA in PDB structures: the bound segments and the methylarginines inside it |
| `supporting/rna_bind_n_seq_v1.6.0.tsv` | methylated protein and ENCODE RNA Bind-n-Seq experiment: its most enriched 5-mers and logos (in vitro) |
| `supporting/dna_motifs_v1.6.0.tsv` | methylated protein and JASPAR 2026 DNA-binding profile |
| `supporting/htr_selex_v1.6.0.tsv` | methylated protein and HTR-SELEX RNA-binding model (Jolma et al. 2020) |
| `supporting/eclip_v1.6.0.tsv` | methylated protein and ENCODE eCLIP experiment: peaks by transcript region and the genes with most peaks (in cells) |
| `supporting/disease_terms_v1.6.0.tsv` | enriched DisGeNET term (clusterProfiler enricher on the DisGeNET terms of DOSE 4.4.0) |
| `supporting/disease_terms_tested_v1.6.0.tsv` | DisGeNET term tested (clusterProfiler enricher on the DisGeNET terms of DOSE 4.4.0), enriched or not |
| `supporting/disease_ontology_v1.6.0.tsv` | Disease Ontology term tested (clusterProfiler enricher on the gene sets of DOSE enrichDO) |
| `supporting/disease_families_v1.6.0.tsv` | disease family: terms and genes |
| `supporting/arginine_background_v1.6.0.tsv` | protein and residue test: the other arginines of the protein, in all and in each class or motif |
| `supporting/density_windows_v1.6.0.tsv` | protein and window width with site density significantly above the permutation null (q <= 0.05) |
| `supporting/condensates_v1.6.0.tsv` | CD-CODE condensate: name, type, members and methylated members |
| `supporting/hpa_locations_v1.6.0.tsv` | protein and Human Protein Atlas location |
| `supporting/other_ptm_sites_v1.6.0.tsv` | other modification on a protein of the release (PTM master table) |
| `supporting/writer_perturbations_v1.6.0.tsv` | site and knockout, knockdown or inhibitor (PRMT1, PRMT5, PRMT7, CARM1): its peptides, change, and the twofold, q < 0.05 and published calls |
| `supporting/functional_experiments_v1.6.0.tsv` | curated experiment on a site or a protein's methylation: perturbation, readout, direction and reference (kept apart from every statistic) |
| `supporting/eclip_rna_compartments_v1.6.0.tsv` | protein with eCLIP data, map and compartment: whether the RNAs it binds gather there, within expression and gene-length deciles (every protein and compartment a map tests) |
| `supporting/methylarginine_site_interaction_records_v1.6.0.tsv` | interaction record at a site: its PDB entry, IntAct interaction, study or model, method, species, paper, the state of the site in the experiment and what the record shows |
| `data_dictionary.tsv` | column of every table, with its codes and their meanings |
| `release.json` | the release: version, counts, reference, inputs, files with MD5 |
| `CHANGES.md` | what changed since the previous version |
| `STUDIES.md` | each study and database behind the sites: where its records come from, how they were filtered and placed, and records traced from the source to the site |
| `MD5SUMS` | checksums (`md5sum -c MD5SUMS`) |

Tables are UTF-8, tab-separated, one header row. Lists are separated by `; `, booleans are `true`/`false`, and a missing value is empty. Coded columns hold codes only; `data_dictionary.tsv` gives their meanings.

### Cite

Cite DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, and the studies of the sites you use (`studies` column, `studies_v1.6.0.tsv`):

| Study or database | Citation | DOI | Records | Sites | Only here |
|---|---|---|---:|---:|---:|
| Larsen 2016 (`Larsen2016`) | Larsen SC et al., Sci Signal 9:rs9 (2016) | [10.1126/scisignal.aaf7329](https://doi.org/10.1126/scisignal.aaf7329) | 17,766 | 7,260 | 1,643 |
| PhosphoSitePlus (PSP) | Hornbeck PV et al., Nucleic Acids Res 43:D512-D520 (2015) | [10.1093/nar/gku1267](https://doi.org/10.1093/nar/gku1267) | 8,375 | 7,443 | 1,229 |
| Maron 2021 (`Maron2021`) | Maron MI et al., iScience 24:102971 (2021) | [10.1016/j.isci.2021.102971](https://doi.org/10.1016/j.isci.2021.102971) | 6,108 | 1,482 | 338 |
| Guo 2014 (`Guo2014`) | Guo A et al., Mol Cell Proteomics 13:372-387 (2014) | [10.1074/mcp.O113.027870](https://doi.org/10.1074/mcp.O113.027870) | 4,729 | 1,569 | 46 |
| Musiani 2019 (`Musiani2019`) | Musiani D et al., Sci Signal 12:eaat8388 (2019) | [10.1126/scisignal.aat8388](https://doi.org/10.1126/scisignal.aat8388) | 4,223 | 1,299 | 226 |
| Fong 2019 (`Fong2019`) | Fong JY et al., Cancer Cell 36:194-209 (2019) | [10.1016/j.ccell.2019.07.003](https://doi.org/10.1016/j.ccell.2019.07.003) | 3,841 | 1,270 | 173 |
| ProMetheusDB | Massignani E et al., Mol Cell Proteomics 21:100243 (2022) | [10.1016/j.mcpro.2022.100243](https://doi.org/10.1016/j.mcpro.2022.100243) | 1,598 | 281 | 77 |
| Li 2021 (`Li2021`) | Li WJ et al., Nat Commun 12:1946 (2021) | [10.1038/s41467-021-21963-1](https://doi.org/10.1038/s41467-021-21963-1) | 2,880 | 319 | 43 |
| Fedoriw 2019 (`Fedoriw2019`) | Fedoriw A et al., Cancer Cell 36:100-114 (2019) | [10.1016/j.ccell.2019.05.014](https://doi.org/10.1016/j.ccell.2019.05.014) | 3,610 | 409 | 29 |
| UniProt | UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025) | [10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010) | 1,016 | 212 | 43 |
| Lim 2020 (`Lim2020`) | Lim Y et al., Proteome Sci 18:6 (2020) | [10.1186/s12953-020-00162-8](https://doi.org/10.1186/s12953-020-00162-8) | 2,001 | 561 | 86 |
| Wei 2020 (`Wei2020`) | Wei M et al., Front Cell Dev Biol 8:545934 (2020) | [10.3389/fcell.2020.545934](https://doi.org/10.3389/fcell.2020.545934) | 578 | 283 | 59 |
| IEDB | Vita R et al., Nucleic Acids Res 47:D339-D343 (2019) | [10.1093/nar/gky1006](https://doi.org/10.1093/nar/gky1006) | 166 | 43 | 20 |
| Hartel 2019 (`Hartel2019`) | Hartel NG et al., Mol Cell Proteomics 18:2149-2164 (2019) | [10.1074/mcp.RA119.001625](https://doi.org/10.1074/mcp.RA119.001625) | 2,806 | 1,447 | 122 |
| Geoghegan 2015 (`Geoghegan2015`) | Geoghegan V et al., Nat Commun 6:6758 (2015) | [10.1038/ncomms7758](https://doi.org/10.1038/ncomms7758) | 4,536 | 2,233 | 203 |
| Sylvestersen 2014 (`Sylvestersen2014`) | Sylvestersen KB et al., Mol Cell Proteomics 13:2072-2088 (2014) | [10.1074/mcp.O113.032748](https://doi.org/10.1074/mcp.O113.032748) | 1,281 | 889 | 29 |
| Lu 2025 (`Lu2025`) | Lu L et al., EMBO Rep 26:5649-5672 (2025) | [10.1038/s44319-025-00599-y](https://doi.org/10.1038/s44319-025-00599-y) | 743 | 155 | 8 |
| Ma 2023 (`Ma2023`) | Ma M et al., J Am Soc Mass Spectrom 34:1692-1700 (2023) | [10.1021/jasms.3c00154](https://doi.org/10.1021/jasms.3c00154) | 403 | 403 | 46 |
| Lu 2024 (`Lu2024`) | Lu L et al., J Proteome Res 23:1014-1027 (2024) | [10.1021/acs.jproteome.3c00724](https://doi.org/10.1021/acs.jproteome.3c00724) | 856 | 398 | 31 |

The site records come from Maron et al. 2021 Table S5 (sites compiled from ten studies), ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of an iPTMnet-format PTM master table (UniProt and IEDB records), and the supplementary tables of Hartel 2019, Geoghegan 2015, Sylvestersen 2014, Lu 2025, Ma 2023, Lu 2024, each staged into one table of records (`scripts/00_stage_study_sites.py`, `inputs/sources/`). Maron et al. 2021 Tables S2 and S3 and the CARM1 knockout data (Shishkova et al. 2017, doi:10.1038/ncomms15571) annotate sites but add none.

### License

The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)); each study keeps its own terms (`licence` in the studies file). The release includes PhosphoSitePlus® (PSP, www.phosphosite.org) records compiled in Maron et al. 2021 Table S5. PhosphoSite, created by Cell Signaling Technology, was licensed under CC BY-NC-SA 3.0 Unported in 2021, when those records were downloaded ([inputs/licences/](../../inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf)). Those terms ask that PSP be cited (Hornbeck et al., Nucleic Acids Res 43:D512-D520, 2015), that bulk downloads list PSP as the source (the `studies` column does), and that web pages showing PSP-derived sites link them to the PSP protein page (www.phosphosite.org/uniprotAccAction?id=<UniProt accession>).

The annotation and interaction sources keep their own terms too (`licence` in `inputs/MANIFEST.tsv`). Five are non-commercial: the DisGeNET terms behind the disease tables (CC BY-NC-SA 4.0), CD-CODE condensates (CC BY-NC 4.0), the pepR-MS RNA crosslink compendium (CC BY-NC), the DNA crosslinks of Trendel et al. 2025 (CC BY-NC 4.0) and the supporting information of Lu et al. 2024 (CC BY-NC 4.0). Maron et al. 2021 and ProMetheusDB supplements are CC BY-NC-ND 4.0: the release redistributes records derived from them, not the tables themselves.

### How it was built

`bash hpc/build_release.sh v1.4.3` at commit `52451320067f754d32574061b8653a2b5f03ec37` of https://github.com/Shechterlab/PTM_Pipeline:

1. checks every input against `inputs/MANIFEST.tsv`;
2. parses the source tables and places each record on the reference sequence (`scripts/02_remap_arg_methyl_to_canonical.py`);
3. stops if the placement of any record contradicts its own evidence (`scripts/qc_remap_evidence_check.py`; the 22 accepted paralog rows of `config/remap_evidence_allowlist.tsv` aside) or any site is not an arginine;
4. runs the analyses on the union: functional classes (`scripts/03`), motifs and arginine enrichment (`08`), disorder and domain context (`06`, `05`, `07`), other PTMs (`09`), site density and clusters (`10`, `28`, `25`), condensates (`12`), Maron et al. 2021 summaries (`16`) and disease enrichment (`14`);
5. encodes the release for the website (`scripts/36`) and computes the interaction evidence of its sites from `inputs/interactions/` (`scripts/37`);
6. writes this folder (`scripts/build_release.py`).

- **Reference sequences:** UniProt release 2026_03 (02-September-2026), the 20,416 reviewed entries of the human reference proteome UP000005640 (`inputs/reference/`). `reference_sequence_md5` tells whether a later UniProt entry still has the same numbering.
- **Context:** MobiDB predicted disorder (export of 2026-03-30), InterPro 108.0 domains and repeats, sequence motifs: the classifiers of `scripts/06`, `scripts/05`, `scripts/07` and `scripts/08`. The intervals were computed on earlier UniProt sequences and moved to these (`inputs/annotation/interval_lift_to_uniprot_2026_03.tsv`).
- **Previous version:** `releases/v1.4.2/` (git tag v1.4.2).

### Annotation layers

Data from other resources, joined to the sites and proteins by accession and position; none adds a site. Each is staged into `inputs/layers/` (the provenance of the interaction records into `inputs/interactions/`; its files listed, with checksums, in `inputs/MANIFEST.tsv`) and registered in `scripts/release_layers.py`, which the build reads for the columns and tests a layer adds (`ADDING_A_SOURCE.md`).

| Layer | Staged by | Adds | License |
|---|---|---|---|
| CARM1 knockout response (Shishkova et al. 2017) | `41_stage_carm1_ko_sites.py` | site column `carm1_ko_call` and `supporting/carm1_ko` (scripts/build_release.py) | CC BY 4.0 |
| Human Protein Atlas 25.1 subcellular locations | `42_stage_hpa_subcellular.py` | protein column `hpa_locations`, `supporting/hpa_locations` and analysis `hpa_location` (scripts/release_tables.py) | CC BY 4.0 |
| Cancer cell lines and tumors | `44b_extract_tcga_mc3.py`, `44_stage_depmap.py` | site columns `tcga_changes`, `tcga_change_tumours`, `depmap_changes`, `depmap_change_lines`; protein columns `depmap_gene_effect_median`, `depmap_lines_screened`, `depmap_dependent_lines`, `depmap_dependency_class`; analyses `depmap_dependency`, `cancer_change` in the statistics table | DepMap Public 24Q4 and PRISM Repurposing 24Q2: CC BY 4.0; TCGA MC3 calls: GDC open-access data (NIH Genomic Data Sharing policy, no license stated) |
| Writer perturbations (PRMT knockouts, knockdowns and inhibitors) | `41_stage_carm1_ko_sites.py`, `45_stage_prmt7_szewczyk2020.py`, `46_stage_hartel2019_prmt1_kd.py`, `51_stage_lu2024_prmt5.py`, `50_stage_li2021_perturbations.py` | site columns `perturbations`, `perturbation_calls`, `perturbation_calls_twofold_q05`, `perturbation_calls_published`; supporting tables `supporting/writer_perturbations` | CC BY 4.0 (Shishkova et al. 2017; Szewczyk et al. 2020; Hartel et al. 2019, from 2020-09-15; Li et al. 2021, its Supplementary Data taken as covered by the article's license); CC BY-NC 4.0 (Lu et al. 2024 supporting information) |
| PaxDb 6.1 protein abundance | `47_stage_paxdb.py` | protein columns `paxdb_abundance_ppm`, `paxdb_percentile` | CC BY 4.0 |
| Reactome 97 pathways | `48_stage_reactome.py` | protein columns `reactome_top_level`, `reactome_n_pathways`; analyses `reactome_pathway`, `reactome_top_level` in the statistics table | CC0 (Reactome data; reactome.org/license) |
| Curated functional experiments | `49_check_functional_curation.py` | supporting tables `supporting/functional_experiments` | CC BY-NC-SA 4.0 (this repository's data license); the facts are from the papers cited |
| Maps of the cell | `52_stage_cell_maps.py` | protein columns `hcm_compartment`, `hcm_domain`, `organelle_ip_compartment`, `chromatin_probability`, `nuclear_body_proteomes`, `eclip_rna_compartments`; analyses `cell_hcm_compartment`, `cell_hcm_domain`, `cell_organelle_ip`, `cell_nuclear_body_proteome`, `cell_chromatin`, `cell_region_hcm`, `cell_region_organelle_ip`, `cell_region_hpa`, `eclip_rna_apexseq`, `eclip_rna_nuclear_body`, `eclip_gene_tsaseq`, `eclip_gene_hic` in the statistics table; supporting tables `supporting/eclip_rna_compartments` | Human Cell Map CC BY 4.0; organelle IP BSD-3; Kustatscher et al. 2014 CC BY; Barutcu et al. 2022 Table S1 (Elsevier): the calls derived from it only; GEO and ENCODE data without restrictions; GENCODE open |
| Provenance of the interaction records | `53_stage_interaction_sources.py`, `37a_stage_interactions.py` | site columns `interaction_classes`, `partners_binding_methylated`, `methylating_enzymes`, `demethylating_enzymes`, `interaction_papers`; supporting tables `supporting/site_interaction_records`; its table and site columns are filled by scripts/build_release.py from supporting/site_interactions | PDB data CC0 1.0 (RCSB Data API and ModelServer); IntAct CC BY 4.0; PSI-MI and PSI-MOD (OLS) CC BY 4.0; PubMed citation data under the NLM Terms and Conditions (courtesy of the U.S. National Library of Medicine) |

CPLM 4.0, listed here in earlier releases for lysine methylation sites drawn on the protein maps, is not a source of this release: no table of it and none of the website data holds CPLM data (CHANGES.md).

Cite them as: CARM1 knockout response (Shishkova et al. 2017): Shishkova et al., Nat Commun 8:15571 (2017); Human Protein Atlas 25.1 subcellular locations: Thul et al., Science 356:eaal3321 (2017); Cancer cell lines and tumors: DepMap Public 24Q4 (doi:10.25452/figshare.plus.27993248; Tsherniak et al. 2017, doi:10.1016/j.cell.2017.06.010; Chronos, Dempster et al. 2021, doi:10.1186/s13059-021-02540-7); PRISM Repurposing 24Q2 (doi:10.6084/m9.figshare.25917643; Corsello et al. 2020, doi:10.1038/s43018-019-0018-6); the results shown here are in whole or part based upon data generated by the TCGA Research Network (https://www.cancer.gov/tcga): MC3 somatic calls (Ellrott et al. 2018, doi:10.1016/j.cels.2018.03.002), the calls of the PanCancer Atlas (Hoadley et al. 2018, doi:10.1016/j.cell.2018.03.022), from the GDC (mc3.v0.2.8.PUBLIC.maf.gz); Writer perturbations (PRMT knockouts, knockdowns and inhibitors): Hartel et al., Mol Cell Proteomics 18:2149-2164 (2019), doi:10.1074/mcp.RA119.001625; Li WJ et al., Nat Commun 12:1946 (2021), doi:10.1038/s41467-021-21963-1; Lu et al., J Proteome Res 23:1014-1027 (2024), doi:10.1021/acs.jproteome.3c00724; Shishkova et al., Nat Commun 8:15571 (2017), doi:10.1038/ncomms15571; Szewczyk et al., Nat Commun 11:2396 (2020), doi:10.1038/s41467-020-16271-z; PaxDb 6.1 protein abundance: PaxDb 6.1, H. sapiens whole-organism integrated dataset (pax-db.org); Huang et al., Mol Cell Proteomics 22:100640 (2023), doi:10.1016/j.mcpro.2023.100640; Reactome 97 pathways: Reactome release 97; Ragueneau et al., Nucleic Acids Res 54:D673-D681 (2026), doi:10.1093/nar/gkaf1223; Curated functional experiments: the paper of each row; Maps of the cell: Go et al., Nature 595:120 (2021), doi:10.1038/s41586-021-03592-2; Hein et al., Cell 188:1137 (2025), doi:10.1016/j.cell.2024.11.028; Kustatscher et al., EMBO J 33:648 (2014), doi:10.1002/embj.201387614; Barutcu et al., Mol Cell 82:1035 (2022), doi:10.1016/j.molcel.2021.12.010; Fazal et al., Cell 178:473 (2019), doi:10.1016/j.cell.2019.05.027; Chen et al., J Cell Biol 217:4025 (2018), doi:10.1083/jcb.201807108; Rao et al., Cell 159:1665 (2014), doi:10.1016/j.cell.2014.11.021; ENCODE Project Consortium, Nature 583:699 (2020), doi:10.1038/s41586-020-2493-4; Provenance of the interaction records: Burley et al., Nucleic Acids Res 51:D488 (2023); Sehnal et al., Nucleic Acids Res 49:W431 (2021); del Toro et al., Nucleic Acids Res 50:D648 (2022); Montecchi-Palazzi et al., Nat Biotechnol 26:864 (2008); McLaughlin et al., Nucleic Acids Res 53:D1318 (2025); Sayers et al., Nucleic Acids Res 53:D20 (2025), doi:10.1093/nar/gkae979.

### Statistics

`statistics_v1.6.0.tsv` holds every test behind the summaries: its analysis and category, the 2x2 counts a, b, c and d, the effect, p, q and the method:

- **Protein and gene sets** (method `enricher`): PANTHER and functional classes, Reactome pathways, CD-CODE condensates, Human Protein Atlas locations, the compartments of the maps of the cell, disease families, Gene Ontology RNA and DNA binding and DepMap dependency classes; the complexes, Disease Ontology terms and DisGeNET terms, tested the same way, are in `supporting/`. Over-representation analysis in R with clusterProfiler::enricher, the DisGeNET terms on DOSE 4.4.0's DisGeNET term-gene table (the test DOSE::enrichDGN ran; `scripts/enrich_sets.R`; DOSE 4.4.0, R 4.5.3, clusterProfiler 4.18.4). For a set holding M of the N members of the background, where n of the N are methylarginine-containing proteins (or their genes) and k of those are in the set, p is the one-sided hypergeometric probability of k or more, `fold_enrichment` is (k/n)/(M/N), and q is the Benjamini-Hochberg adjusted p over the sets of the analysis; a = k, b = n - k, c = M - k and d = N - M - n + k. A set is tested when it holds a methylarginine-containing protein and its size is within the analysis's range. Backgrounds: the reviewed human proteome (classes, Gene Ontology, complexes), the proteins of any human CD-CODE condensate (condensates), the proteins the Human Protein Atlas profiled (locations), the reviewed human genes with a DisGeNET term (disease families and DisGeNET terms) and those with a Disease Ontology term (Disease Ontology), the reviewed proteins whose gene DepMap 24Q4 screened (dependency classes), and the proteins each map of the cell places (analyses cell_*). The analyses cell_region_* group the compartments of the Human Cell Map, organelle IP and the Human Protein Atlas into three regions of the cell, each tested as one set: the nucleus, the cytoplasm outside the organelles, and the organelles and membranes with the nuclear envelope (`scripts/cell_maps.py`, REGION_OF).
- **Proteins, within abundance deciles** (method `mantel_haenszel`, the background ending "each PaxDb abundance decile a stratum"): mass spectrometry finds abundant proteins most, so a set rich in abundant proteins can look enriched among the methylated proteins whatever methylation does. Each protein-level over-representation test above (PANTHER, functional classes, Reactome, Human Protein Atlas locations, the maps of the cell, condensates, DepMap dependency, Gene Ontology binding, disease families, and the proteins found binding RNA or DNA in cells and the eCLIP'd proteins below) is run again with its background split into the deciles of PaxDb 6.1 abundance, plus one stratum for the proteins without a PaxDb value: the Mantel-Haenszel common odds ratio with its Robins-Breslow-Greenland 95% CI (statsmodels) and the continuity-corrected two-sided Cochran-Mantel-Haenszel test, its |O - E| - 0.5 floored at zero (`common.cmh_statistic`), q Benjamini-Hochberg within the analysis. The sets and the pooled a, b, c, d are those of the enricher row beside it; genes (disease families) take the highest PaxDb value of their reviewed entries. This is a sensitivity test, reported beside the test it repeats: it asks whether an enrichment holds among proteins of similar abundance, and removes any real association of methylation with abundance too (`scripts/abundance_strata.py`).
- **Arginines** (methods `fisher` and `mantel_haenszel`): motifs, disorder, domain context, nearest domain class, structural context, the interaction tests and the changes in tumors and cell lines (cancer_change) compare methylated arginines with the other arginines of the same proteins: two-sided Fisher exact test (scipy) with the sample odds ratio after adding 0.5 to each cell, pooled over the arginines of the methylated proteins, and the Mantel-Haenszel common odds ratio with each protein a stratum (statsmodels), which compares each protein's methylated arginines with its own other arginines; q is Benjamini-Hochberg within the analysis and method.
- **Proteins found binding RNA or DNA in cells** (`interaction_statistics`, analysis binding_protein, methods `enricher` and `mantel_haenszel`): methylarginine-containing proteins against the other reviewed human proteins among the proteins each experiment measured (the reviewed proteome where it lists none), tested as Gene Ontology binding is: clusterProfiler enricher, one job per experiment since each has its own background, the one-sided hypergeometric p and fold enrichment, q Benjamini-Hochberg over the experiments; and again within PaxDb abundance deciles (above), since pepR-MS and the DNA crosslinks list no proteins measured and abundance can matter beyond detection.
- **Where the RNAs bound by proteins with eCLIP data gather** (analyses eclip_*, scripts/cell_maps.py): for each protein with eCLIP data and each map of RNA or DNA (APEX-seq, nuclear-body APEX2-seq, TSA-seq, Hi-C subcompartments), its bound RNAs against the other RNAs (genes) the map places, in each compartment: the Mantel-Haenszel odds ratio within deciles of expression in the map's cells crossed with deciles of gene length (GENCODE v19; the genes it lacks are one more length stratum), q Benjamini-Hochberg over every protein and compartment of the map (`supporting/eclip_rna_compartments`; its RNAs gather there at q <= 0.05 and odds ratio above 1). The methylated proteins are then tested against the other eCLIP'd proteins the map tests, as a protein set (methods `enricher` and `mantel_haenszel`, above), and the size of every protein's odds ratio, called or not, is compared too (method `mann_whitney`: two-sided Mann-Whitney U, a and c the methylated and the other proteins compared, q Benjamini-Hochberg within the map).

### Columns

The same definitions, with every code and its meaning, are in `data_dictionary.tsv`.

#### `human_methylarginine_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | UniProt accession and arginine position, e.g. P35637:R216. |
| `uniprot_accession` | string | reviewed UniProtKB accession of the canonical sequence. |
| `uniprot_entry_name` | string | UniProtKB entry name. |
| `gene` | list&lt;string&gt; | primary gene name of the UniProt entry (several for a few entries, e.g. histone genes), empty when UniProt gives none. |
| `protein_name` | string | recommended protein name of the UniProt entry. |
| `position` | integer | position of the arginine in the reference sequence (1-based, initiator methionine = 1). |
| `residue` | code | residue at that position in the reference sequence: always R. Codes: `R`. |
| `window_15` | string | residues -7 to +7 around the site, _ beyond a protein end. |
| `window_31` | string | residues -15 to +15 around the site, _ beyond a protein end. |
| `protein_length` | integer | length of the reference sequence. |
| `reference_sequence_md5` | string | MD5 of the reference sequence, uppercase hex as UniProt gives it (a later UniProt entry with the same MD5 has the same numbering). |
| `methyl_states_reported` | list&lt;code&gt; | methylation degrees given by the source records, in code order, empty where they give none. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `maron2021_s2_states` | list&lt;code&gt; | degrees called at the site in Maron et al. 2021 Table S2. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `uniprot_states` | list&lt;code&gt; | degrees UniProt release 2026_03 annotates at the site (MOD_RES) with experimental evidence: ECO:0000269, or ECO:0007744 (large-scale proteomics). ADMA and SDMA calls in UniProt mostly follow the antibody of the experiment. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `uniprot_states_inferred` | list&lt;code&gt; | degrees UniProt annotates at the site only by similarity to another species (ECO:0000250) or without experimental evidence; not counted in methyl_state. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `states_by_study_design` | list&lt;code&gt; | the one degree a study looked for: Rme1 for Larsen 2016 (monomethylarginine antibodies; Table S1 lists monomethylarginine sites, and the deposited search, PRIDE PXD003700, has monomethylation only), Rme2 for Wei 2020 (ADMA antibody, PTMScan #13474; only dimethylation searched). Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `methyl_state` | code | the degrees reported at the site by any of methyl_states_reported, maron2021_s2_states, uniprot_states and states_by_study_design; empty when none gives one. Mass spectrometry distinguishes mono- from dimethylation; asymmetric and symmetric dimethylation are mostly assigned from the antibody used. Codes: `Rme1`, `Rme2`, `Rme1_Rme2`. |
| `methyl_state_evidence` | list&lt;code&gt; | which of these give a degree, in code order. Codes: `source_records`, `maron2021_s2`, `uniprot`, `study_design`. |
| `maron2021_s3_antibodies` | list&lt;code&gt; | antibodies whose immunoprecipitation recovered the site in Maron et al. 2021 Table S3 (supporting/maron2021_detection). Codes: `MMA`, `ADMA`, `SDMA`. |
| `dimethyl_type_neutral_loss` | list&lt;code&gt; | the dimethylarginine type Hartel et al. 2019 called at the site from neutral losses in the fragment spectra, independent of the antibody used (both when both were called). Shown beside maron2021_s3_antibodies, which says which antibodies recovered the site; the two can differ. Codes: `ADMA`, `SDMA`. |
| `maron2021_GSK591_call` | code | response to GSK591 in Maron et al. 2021 Table S3, combined over antibodies (supporting/maron2021_inhibitors: a peptide changes at adjusted p <= 0.05): peptide-level evidence, as a peptide carrying several methylarginines gives each the same values (the supporting table separates single-site peptides); empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `maron2021_MS023_call` | code | response to MS023, as for GSK591. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `carm1_ko_call` | code | response to CARM1 knockout (Shishkova et al. 2017), combined over cell lines (supporting/carm1_ko: twofold at BH q < 0.05): peptide-level evidence; empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `motifs` | list&lt;code&gt; | sequence motifs around the site (scripts/08), in code order. Codes: `RG`, `RGG`, `GRG`, `GAR`, `RXR`, `CARM1_proline_rich`, `CARM1_hydrophobic_proline`, `DR`, `RD`, `ER`, `RE`. |
| `disorder` | code | MobiDB predicted disorder at the site (scripts/06). Codes: `disordered`, `disorder_boundary`, `ordered`, `no_disorder_annotation`. |
| `disorder_edge_distance` | integer | residues to the nearest edge of a disordered region (a protein end is not an edge), empty with no annotation or when the region spans the whole protein. |
| `domain` | code | InterPro 108.0 domain or repeat context (scripts/05). Codes: `in_domain`, `boundary`, `inter_domain_linker`, `distal`, `no_domain_annotation`. |
| `nearest_domain_interpro` | string | InterPro accession of the nearest domain or repeat (the innermost when inside several). |
| `nearest_domain_name` | string | its InterPro name. |
| `nearest_domain_start` | integer | its first residue. |
| `nearest_domain_end` | integer | its last residue. |
| `nearest_domain_class` | string | its class in config/domain_class_ontology.json, as scripts/07 assigns it, No domain annotation when none. |
| `domain_side` | code | where the site lies relative to the nearest domain, empty when none. Codes: `inside`, `N_terminal_to_domain`, `C_terminal_to_domain`. |
| `domain_edge_distance` | integer | residues to the nearest domain edge (0 inside a domain), empty when none. |
| `structure_context` | code | disorder and domains in one class, from inside a folded domain to deep in a disordered region (common.structure_context, from disorder, domain and domain_edge_distance). Codes: `domain_folded`, `domain_disordered`, `ordered_outside_domain`, `idr_near_domain`, `idr_mid_domain`, `idr_far_from_domain`, `no_disorder_annotation`. |
| `interaction_status` | code | interaction evidence at the site (supporting/site_interactions), empty for the sites kept pending. Codes: `experimental`, `predicted_only`, `none`. |
| `interaction_evidence` | list&lt;code&gt; | the evidence types at the site, in code order. Codes: see `data_dictionary.tsv`. |
| `bound_nucleic_acid` | list&lt;code&gt; | RNA or DNA in direct contact with the site in PDB structures, within 8 Å of its side chain in a complex without direct contact, or in contact via a residue at i ± 1-5 in a complex with neither, in code order (every entry: supporting/nucleic_acid_contacts and nucleic_acid_near_sites). Codes: `rna_at`, `dna_at`, `rna_proximal`, `dna_proximal`, `rna_near`, `dna_near`. |
| `bound_nucleic_acid_types` | list&lt;code&gt; | the types of those RNA and DNA chains, in code order. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `bound_sequence` | string | the sequence window with the most support (direct contact, then proximal, then via a flanking residue; then the most PDB entries): 8 nucleotides on each side of the center, the contacted nucleotides in brackets, - beyond the chain's ends. It is the sequence deposited in the structure, often chosen by its authors, and records a contact, not a sequence preference. |
| `bound_sequence_offset` | integer | 0 for direct contact or proximity of the arginine side chain; else k, the offset of the nearest residue i+k (k = ±1 to ±5) in direct contact (within 4 Å) with bound_sequence. |
| `bound_sequence_distance_a` | number | minimum heavy-atom distance from the arginine side chain to bound_sequence, in Å (direct and proximal), empty for a flanking-residue contact. |
| `bound_sequence_type` | code | the type of its chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `bound_sequence_molecule` | string | its molecule name in the PDB entries (the most common one). |
| `bound_sequence_entries` | integer | number of PDB entries containing it. |
| `tcga_changes` | list&lt;string&gt; | missense and stop changes of the arginine in TCGA tumors (MC3 public calls from the GDC, 33 projects, 10,294 tumor samples; FILTER PASS, and in OV and LAML also wga alone, as MC3 used them), most tumors first, * a stop. A change is read on a transcript UniProt assigns to the canonical isoform where the variant falls in one, else on the MAF's transcript, and counts only where the canonical sequence has this arginine at that position (inputs/layers/README.md). |
| `tcga_change_tumours` | list&lt;integer&gt; | for each of these changes, the tumor samples with it, same order. |
| `depmap_changes` | list&lt;string&gt; | missense and stop changes of the arginine in DepMap 24Q4 cell lines (OmicsSomaticMutations, on the canonical UniProt isoform), most lines first. |
| `depmap_change_lines` | list&lt;integer&gt; | for each of these changes, the cell lines with it, same order. |
| `perturbations` | list&lt;code&gt; | the knockouts, knockdowns and inhibitors that quantified a peptide carrying the site (or, for one that lists only its changed sites, listed it), in code order (supporting/writer_perturbations). Codes: `carm1_ko_mcf7`, `carm1_ko_mdamb231`, `prmt7_ko_hct116`, `prmt1_kd_hek293t`, `prmt5i_epz015666_hela`, `prmt7_kd_hek293`, `prmt7i_sgc3027_hek293`, `carm1_kd_hek293`, `prmt5_kd_hek293`. |
| `perturbation_calls` | list&lt;code&gt; | for each of them, the call (twofold at q < 0.05, or the published call where its code says so): peptide-level evidence, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `perturbation_calls_twofold_q05` | list&lt;code&gt; | for each of them, the call at |log2| >= 1 and q < 0.05, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `perturbation_calls_published` | list&lt;code&gt; | for each of them, the published call, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `interaction_classes` | list&lt;code&gt; | the classes of the interaction records at the site (supporting/site_interaction_records), in code order. Codes: see `data_dictionary.tsv`. |
| `partners_binding_methylated` | list&lt;string&gt; | partners recorded binding the site methylated, by gene (the species in brackets when not human; the name when there is no gene), sorted. |
| `methylating_enzymes` | list&lt;string&gt; | enzymes IntAct records methylating the site (methylation reaction), the same way. |
| `demethylating_enzymes` | list&lt;string&gt; | enzymes IntAct records demethylating the site (demethylation reaction), the same way. |
| `interaction_papers` | integer | distinct papers (PubMed or DOI) behind the experimental records at the site (predicted models left out). |
| `n_studies` | integer | number of studies or databases reporting the site. |
| `studies` | list&lt;code&gt; | the studies or databases reporting the site, in the order of the studies file, PhosphoSitePlus is PSP, as its terms ask. Codes: see `data_dictionary.tsv`. |
| `study_dois` | list&lt;string&gt; | their DOIs, in the same order. |
| `source_tables` | list&lt;code&gt; | the source tables of its records, in code order. Codes: see `data_dictionary.tsv`. |
| `n_records` | integer | number of source records placed on the site. |
| `record_ids` | list&lt;integer&gt; | the records behind the site (record_id in supporting/source_records), ascending. |
| `localised_by` | list&lt;code&gt; | how its records localize the methyl group to this arginine, in code order: the kinds of evidence behind the site, which do not share one error rate. Codes: `localisation_score`, `one_methylatable_residue`, `every_arginine_dimethylated`, `curated_experimental`, `curated_database`. |
| `curated_only` | boolean | true when every record placing the site is a curated annotation (UniProt with experimental evidence, PhosphoSitePlus as compiled in 2021), with no localization score or sequence behind it. |
| `reported_positions` | list&lt;integer&gt; | positions the records report, ascending, only when one differs from position. |
| `flags` | list&lt;code&gt; | cautions, in code order. Codes: `one_sided_window`, `paralog_window`, `near_match_only`, `arginine_from_position`, `paralog_assigned`, `shared_sequence`, `single_study`, `xenograft_only`. |
| `shared_with_paralogs` | list&lt;string&gt; | for a site flagged paralog_assigned: the other reviewed entries whose sequence holds the same peptide or window at an arginine, sorted. |
| `same_window_accessions` | list&lt;string&gt; | every other reviewed entry holding the 31 residues of window_31, sorted. |
| `same_window_genes` | list&lt;string&gt; | the first primary gene of each (the accession when UniProt gives no gene), same order. |
| `in_previous_version` | boolean | whether the previous version had this site. |
| `change_since_previous_version` | code | kept or added (removed sites are in supporting/site_changes). Codes: `kept`, `added`, `removed`. |
| `change_reasons` | list&lt;code&gt; | why the site was added, or renumbered (CHANGES.md). Codes: see `data_dictionary.tsv`. |

#### `human_methylarginine_proteins_v1.6.0.tsv`

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
| `studies` | list&lt;code&gt; | studies or databases reporting its sites, in the order of the studies file. Codes: see `data_dictionary.tsv`. |
| `panther_group` | string | group of its PANTHER 19.0 top-level protein class (config/panther_class_groups.json, scripts/03b), or No PANTHER protein class. |
| `panther_class` | string | its PANTHER top-level protein class, empty when PANTHER gives none. |
| `panther_classes` | list&lt;string&gt; | its most specific PANTHER protein classes. |
| `panther_family` | string | its PANTHER family (PTHR id). |
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
| `disease_family_terms` | list&lt;integer&gt; | for each of these families, the number of its DisGeNET terms annotated to the gene. |
| `als_ftd_family_gene` | boolean | the gene is in the ALS / motor neuron / FTD family. |
| `disgenet_background` | boolean | the gene has at least one DisGeNET term: it is in the background of the disease family tests. |
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
| `depmap_gene_effect_median` | number | median CRISPR gene effect (Chronos) of the protein's gene across the DepMap 24Q4 lines that screened it: 0 no effect, -1 the median of commonly essential genes; empty when DepMap did not screen the gene (an isoform takes its entry's gene). |
| `depmap_lines_screened` | integer | lines with a gene effect for the gene. |
| `depmap_dependent_lines` | integer | lines whose dependency probability for the gene is at least 0.5. |
| `depmap_dependency_class` | code | the gene's dependency class in DepMap 24Q4, empty when not screened. Codes: `common_essential`, `most_lines`, `some_lines`, `few_or_none`. |
| `paxdb_abundance_ppm` | number | protein abundance in parts per million of the proteome, PaxDb 6.1 H. sapiens whole-organism integrated dataset (a weighted average of 7 whole-organism datasets, PeptideAtlas 2011 weighing 0.75), through PaxDb's UniProt cross-reference to one STRING v12.0 protein (else UniProt's STRING cross-reference); empty when PaxDb has no value for the entry (an isoform takes its entry's value). For a detection covariate, not a measurement in the cells the sites come from. |
| `paxdb_percentile` | number | percentile of that abundance among the 18,805 reviewed entries with a value (ties averaged; 100 the most abundant). |
| `reactome_top_level` | list&lt;string&gt; | the Reactome 97 top-level pathways the protein is in, at any level below them, sorted (an isoform takes its entry's); empty when Reactome gives the entry no pathway. |
| `reactome_n_pathways` | integer | the entry's lowest-level Reactome 97 pathways (UniProt2Reactome); empty when Reactome gives the entry none. |
| `hcm_compartment` | string | the Human Cell Map compartment of the entry (Go et al. 2021, HEK293 BioID, NMF localization); empty when the map does not localize it. |
| `hcm_domain` | string | its Human Cell Map SAFE domain, the finer neighbourhood of proteins that label one another; empty when SAFE gives it none. |
| `organelle_ip_compartment` | string | its organelle IP cluster (Hein et al. 2025, HEK293T), or unclassified; empty when organelle IP did not detect it. |
| `chromatin_probability` | number | its interphase chromatin probability (Kustatscher et al. 2014): how reliably it comes down with chromatin; empty when not scored. |
| `nuclear_body_proteomes` | list&lt;string&gt; | the nuclear speckles, lamina and PML bodies whose APEX2 proximity proteome holds it (Barutcu et al. 2022, SAINT AvgP >= 0.9; speckles with two of the three markers); none when detected in none; empty when not detected. |
| `eclip_rna_compartments` | list&lt;string&gt; | where the RNAs its eCLIP experiments bind gather, as map: compartment (q <= 0.05 and odds ratio above 1 within expression and gene-length deciles; supporting/eclip_rna_compartments); none when a map tests it and they gather nowhere; empty when no map tests it (no eCLIP data, or fewer than 20 bound RNAs placed). |

#### `human_methylarginine_lower_confidence_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | UniProt accession and arginine position, e.g. P35637:R216. |
| `uniprot_accession` | string | reviewed UniProtKB accession of the canonical sequence. |
| `uniprot_entry_name` | string | UniProtKB entry name. |
| `gene` | list&lt;string&gt; | primary gene name of the UniProt entry (several for a few entries, e.g. histone genes), empty when UniProt gives none. |
| `protein_name` | string | recommended protein name of the UniProt entry. |
| `position` | integer | position of the arginine in the reference sequence (1-based, initiator methionine = 1). |
| `residue` | code | residue at that position in the reference sequence: always R. Codes: `R`. |
| `window_15` | string | residues -7 to +7 around the site, _ beyond a protein end. |
| `window_31` | string | residues -15 to +15 around the site, _ beyond a protein end. |
| `protein_length` | integer | length of the reference sequence. |
| `reference_sequence_md5` | string | MD5 of the reference sequence, uppercase hex as UniProt gives it (a later UniProt entry with the same MD5 has the same numbering). |
| `methyl_states_reported` | list&lt;code&gt; | methylation degrees given by the source records, in code order, empty where they give none. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `maron2021_s2_states` | list&lt;code&gt; | degrees called at the site in Maron et al. 2021 Table S2. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `uniprot_states` | list&lt;code&gt; | degrees UniProt release 2026_03 annotates at the site (MOD_RES) with experimental evidence: ECO:0000269, or ECO:0007744 (large-scale proteomics). ADMA and SDMA calls in UniProt mostly follow the antibody of the experiment. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `uniprot_states_inferred` | list&lt;code&gt; | degrees UniProt annotates at the site only by similarity to another species (ECO:0000250) or without experimental evidence; not counted in methyl_state. Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `states_by_study_design` | list&lt;code&gt; | the one degree a study looked for: Rme1 for Larsen 2016 (monomethylarginine antibodies; Table S1 lists monomethylarginine sites, and the deposited search, PRIDE PXD003700, has monomethylation only), Rme2 for Wei 2020 (ADMA antibody, PTMScan #13474; only dimethylation searched). Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `methyl_state` | code | the degrees reported at the site by any of methyl_states_reported, maron2021_s2_states, uniprot_states and states_by_study_design; empty when none gives one. Mass spectrometry distinguishes mono- from dimethylation; asymmetric and symmetric dimethylation are mostly assigned from the antibody used. Codes: `Rme1`, `Rme2`, `Rme1_Rme2`. |
| `methyl_state_evidence` | list&lt;code&gt; | which of these give a degree, in code order. Codes: `source_records`, `maron2021_s2`, `uniprot`, `study_design`. |
| `maron2021_s3_antibodies` | list&lt;code&gt; | antibodies whose immunoprecipitation recovered the site in Maron et al. 2021 Table S3 (supporting/maron2021_detection). Codes: `MMA`, `ADMA`, `SDMA`. |
| `dimethyl_type_neutral_loss` | list&lt;code&gt; | the dimethylarginine type Hartel et al. 2019 called at the site from neutral losses in the fragment spectra, independent of the antibody used (both when both were called). Shown beside maron2021_s3_antibodies, which says which antibodies recovered the site; the two can differ. Codes: `ADMA`, `SDMA`. |
| `maron2021_GSK591_call` | code | response to GSK591 in Maron et al. 2021 Table S3, combined over antibodies (supporting/maron2021_inhibitors: a peptide changes at adjusted p <= 0.05): peptide-level evidence, as a peptide carrying several methylarginines gives each the same values (the supporting table separates single-site peptides); empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `maron2021_MS023_call` | code | response to MS023, as for GSK591. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `carm1_ko_call` | code | response to CARM1 knockout (Shishkova et al. 2017), combined over cell lines (supporting/carm1_ko: twofold at BH q < 0.05): peptide-level evidence; empty when not measured. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `motifs` | list&lt;code&gt; | sequence motifs around the site (scripts/08), in code order. Codes: `RG`, `RGG`, `GRG`, `GAR`, `RXR`, `CARM1_proline_rich`, `CARM1_hydrophobic_proline`, `DR`, `RD`, `ER`, `RE`. |
| `disorder` | code | MobiDB predicted disorder at the site (scripts/06). Codes: `disordered`, `disorder_boundary`, `ordered`, `no_disorder_annotation`. |
| `disorder_edge_distance` | integer | residues to the nearest edge of a disordered region (a protein end is not an edge), empty with no annotation or when the region spans the whole protein. |
| `domain` | code | InterPro 108.0 domain or repeat context (scripts/05). Codes: `in_domain`, `boundary`, `inter_domain_linker`, `distal`, `no_domain_annotation`. |
| `nearest_domain_interpro` | string | InterPro accession of the nearest domain or repeat (the innermost when inside several). |
| `nearest_domain_name` | string | its InterPro name. |
| `nearest_domain_start` | integer | its first residue. |
| `nearest_domain_end` | integer | its last residue. |
| `nearest_domain_class` | string | its class in config/domain_class_ontology.json, as scripts/07 assigns it, No domain annotation when none. |
| `domain_side` | code | where the site lies relative to the nearest domain, empty when none. Codes: `inside`, `N_terminal_to_domain`, `C_terminal_to_domain`. |
| `domain_edge_distance` | integer | residues to the nearest domain edge (0 inside a domain), empty when none. |
| `structure_context` | code | disorder and domains in one class, from inside a folded domain to deep in a disordered region (common.structure_context, from disorder, domain and domain_edge_distance). Codes: `domain_folded`, `domain_disordered`, `ordered_outside_domain`, `idr_near_domain`, `idr_mid_domain`, `idr_far_from_domain`, `no_disorder_annotation`. |
| `tcga_changes` | list&lt;string&gt; | missense and stop changes of the arginine in TCGA tumors (MC3 public calls from the GDC, 33 projects, 10,294 tumor samples; FILTER PASS, and in OV and LAML also wga alone, as MC3 used them), most tumors first, * a stop. A change is read on a transcript UniProt assigns to the canonical isoform where the variant falls in one, else on the MAF's transcript, and counts only where the canonical sequence has this arginine at that position (inputs/layers/README.md). |
| `tcga_change_tumours` | list&lt;integer&gt; | for each of these changes, the tumor samples with it, same order. |
| `depmap_changes` | list&lt;string&gt; | missense and stop changes of the arginine in DepMap 24Q4 cell lines (OmicsSomaticMutations, on the canonical UniProt isoform), most lines first. |
| `depmap_change_lines` | list&lt;integer&gt; | for each of these changes, the cell lines with it, same order. |
| `perturbations` | list&lt;code&gt; | the knockouts, knockdowns and inhibitors that quantified a peptide carrying the site (or, for one that lists only its changed sites, listed it), in code order (supporting/writer_perturbations). Codes: `carm1_ko_mcf7`, `carm1_ko_mdamb231`, `prmt7_ko_hct116`, `prmt1_kd_hek293t`, `prmt5i_epz015666_hela`, `prmt7_kd_hek293`, `prmt7i_sgc3027_hek293`, `carm1_kd_hek293`, `prmt5_kd_hek293`. |
| `perturbation_calls` | list&lt;code&gt; | for each of them, the call (twofold at q < 0.05, or the published call where its code says so): peptide-level evidence, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `perturbation_calls_twofold_q05` | list&lt;code&gt; | for each of them, the call at |log2| >= 1 and q < 0.05, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `perturbation_calls_published` | list&lt;code&gt; | for each of them, the published call, same order. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `n_studies` | integer | number of studies or databases reporting the site. |
| `studies` | list&lt;code&gt; | the studies or databases reporting the site, in the order of the studies file, PhosphoSitePlus is PSP, as its terms ask. Codes: see `data_dictionary.tsv`. |
| `study_dois` | list&lt;string&gt; | their DOIs, in the same order. |
| `source_tables` | list&lt;code&gt; | the source tables of its records, in code order. Codes: see `data_dictionary.tsv`. |
| `n_records` | integer | number of source records placed on the site. |
| `record_ids` | list&lt;integer&gt; | the records behind the site (record_id in supporting/source_records), ascending. |
| `localised_by` | list&lt;code&gt; | how its records localize the methyl group to this arginine, in code order: the kinds of evidence behind the site, which do not share one error rate. Codes: `localisation_score`, `one_methylatable_residue`, `every_arginine_dimethylated`, `curated_experimental`, `curated_database`. |
| `curated_only` | boolean | true when every record placing the site is a curated annotation (UniProt with experimental evidence, PhosphoSitePlus as compiled in 2021), with no localization score or sequence behind it. |
| `reported_positions` | list&lt;integer&gt; | positions the records report, ascending, only when one differs from position. |
| `flags` | list&lt;code&gt; | cautions, in code order. Codes: `one_sided_window`, `paralog_window`, `near_match_only`, `arginine_from_position`, `paralog_assigned`, `shared_sequence`, `single_study`, `xenograft_only`. |
| `shared_with_paralogs` | list&lt;string&gt; | for a site flagged paralog_assigned: the other reviewed entries whose sequence holds the same peptide or window at an arginine, sorted. |
| `same_window_accessions` | list&lt;string&gt; | every other reviewed entry holding the 31 residues of window_31, sorted. |
| `same_window_genes` | list&lt;string&gt; | the first primary gene of each (the accession when UniProt gives no gene), same order. |
| `in_analyses_with_lower_confidence` | boolean | counted in the analyses with lower-confidence sites (statistics_with_lower_confidence): a record covers this arginine with its own evidence (evidence_at_site pass) and leaves the methyl group on an arginine (methyl_on_arginine not fail). |
| `test_evidence_at_site` | code | site test evidence_at_site, for the record that passes most tests. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_named_entry` | code | site test named_entry, same record. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_unique_in_proteome` | code | site test unique_in_proteome, same record. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_methyl_on_arginine` | code | site test methyl_on_arginine, same record. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_localised` | code | site test localized, same record. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `tests_failed` | list&lt;code&gt; | the tests that record fails, in code order. Codes: `evidence_at_site`, `named_entry`, `unique_in_proteome`, `methyl_on_arginine`, `localised`. |
| `reasons` | list&lt;code&gt; | why its records add no site, in code order. Codes: `below_threshold`, `several_residues`, `search_disagrees`, `no_localising_evidence`, `window_fits_several`, `several_sites`, `other_protein`, `not_in_reference`. |
| `in_earlier_versions` | list&lt;string&gt; | the earlier releases that had this site (v1.0.0, v1.1.0). |

#### `studies_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `study_id` | code | study or database. Codes: see `data_dictionary.tsv`. |
| `label` | string | name to show. |
| `source_table` | list&lt;code&gt; | the source tables its records come from. Codes: see `data_dictionary.tsv`. |
| `citation` | string | citation, checked against Crossref. |
| `doi` | string | DOI. |
| `licence` | string | terms of the records as this release received them. |
| `n_records` | integer | source records. |
| `n_records_placed` | integer | of these, placed on a site. |
| `n_sites` | integer | sites (status site) it reports. |
| `n_sites_only_here` | integer | sites no other study or database reports. |
| `n_sites_in_previous_version` | integer | of its sites, those the previous version had. |
| `n_records_several_sites` | integer | records not placed because their evidence names several sites. |
| `n_records_other_protein` | integer | records not placed because their evidence names one site, in a protein the record does not name. |
| `n_records_not_localised` | integer | records not placed because the source gives no localization of the methyl group (placement superseded_by_study_table or unresolved_not_localised). |
| `n_records_no_site` | integer | records not placed for another reason (placement_code in supporting/source_records). |
| `note` | string |  |
| `quality` | string | what a record holds, how the study localized the methyl group, and what was kept. |

#### `statistics_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `analysis` | code | what was tested. Codes: see `data_dictionary.tsv`. |
| `category` | string | the class, motif, family, condensate (CD-CODE id) or location tested. |
| `stratum` | string | all, or the subset the test is restricted to. |
| `a` | integer | targets in the category. |
| `b` | integer | targets not in the category. |
| `c` | integer | background in the category. |
| `d` | integer | background not in the category. |
| `odds_ratio` | number | odds ratio, for the 2x2 tests (methods fisher, mantel_haenszel); empty for method enricher. |
| `log2_odds_ratio` | number | its log2. |
| `fold_enrichment` | number | fold enrichment (k/n over M/N) as clusterProfiler reports it, for method enricher; empty otherwise. |
| `log2_fold_enrichment` | number | its log2. |
| `p` | number | p value. |
| `q` | number | Benjamini-Hochberg adjusted p value within the analysis and method, over all its strata and categories. |
| `method` | code | the test. Codes: `enricher`, `fisher`, `mantel_haenszel`, `mann_whitney`. |
| `background` | string | what the targets are compared with. |
| `n_strata` | integer | strata of a stratified test, empty otherwise. |
| `odds_ratio_ci_low` | number | lower bound of the 95% confidence interval of the odds ratio (Robins-Breslow-Greenland, for the mantel_haenszel rows; empty where the odds ratio is infinite, 0 or not estimable), empty otherwise. |
| `odds_ratio_ci_high` | number | upper bound, as for odds_ratio_ci_low. |

#### `statistics_with_lower_confidence_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `analysis` | code | what was tested. Codes: see `data_dictionary.tsv`. |
| `category` | string | the class, motif, family, condensate (CD-CODE id) or location tested. |
| `stratum` | string | all, or the subset the test is restricted to. |
| `a` | integer | targets in the category. |
| `b` | integer | targets not in the category. |
| `c` | integer | background in the category. |
| `d` | integer | background not in the category. |
| `odds_ratio` | number | odds ratio, for the 2x2 tests (methods fisher, mantel_haenszel); empty for method enricher. |
| `log2_odds_ratio` | number | its log2. |
| `fold_enrichment` | number | fold enrichment (k/n over M/N) as clusterProfiler reports it, for method enricher; empty otherwise. |
| `log2_fold_enrichment` | number | its log2. |
| `p` | number | p value. |
| `q` | number | Benjamini-Hochberg adjusted p value within the analysis and method, over all its strata and categories. |
| `method` | code | the test. Codes: `enricher`, `fisher`, `mantel_haenszel`, `mann_whitney`. |
| `background` | string | what the targets are compared with. |
| `n_strata` | integer | strata of a stratified test, empty otherwise. |
| `odds_ratio_ci_low` | number | lower bound of the 95% confidence interval of the odds ratio (Robins-Breslow-Greenland, for the mantel_haenszel rows; empty where the odds ratio is infinite, 0 or not estimable), empty otherwise. |
| `odds_ratio_ci_high` | number | upper bound, as for odds_ratio_ci_low. |

#### `supporting/source_records_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `record_id` | integer | row number of the record in the source tables as the pipeline reads them, the same in every release. |
| `source_table` | code | the source table the record comes from. Codes: see `data_dictionary.tsv`. |
| `study` | code | the study or database that reported it. Codes: see `data_dictionary.tsv`. |
| `study_doi` | string | its DOI. |
| `reference_text` | string | the reference as the source table gives it. |
| `pmid` | string | PubMed ID, where the source gives one. |
| `reported_accession` | list&lt;string&gt; | accession as reported (a few records give two). |
| `reported_gene` | string | gene as reported. |
| `reported_position` | string | position as reported. |
| `reported_methyl_state` | string | methylation state as reported, in the codes of the site table (ProMetheusDB me = Rme1, di = Rme2). |
| `reported_enzyme` | string | methyltransferase as reported (PTM master table only). |
| `reported_sample` | list&lt;string&gt; | cells or tissue the study measured (records read from a study table). |
| `reported_enrichment` | list&lt;string&gt; | how the study enriched the peptides, e.g. ADMA IAP (antibody), SCX fractionation, SMN Tudor domain (records read from a study table). |
| `reported_detail` | string | what else the study gives for the record: the neutral-loss call, the localization probability or score, the source rows merged into the record and their tables, the PRMT1 knockdown result (records read from a study table). |
| `evidence_kind` | code | kind of sequence evidence. Codes: `window_13`, `window_15`, `window_31`, `peptide`, `marked_peptide`, `position_only`. |
| `evidence_sequence` | list&lt;string&gt; | the window, peptide or marked peptide (a few records give several), as the pipeline reads it: modifications other than methylation, e.g. (ox), are left out, and source_correction says when it was corrected. |
| `evidence_window` | list&lt;string&gt; | the sequence window a marked-peptide record also carries (ProMetheusDB SeqWindow), which the placement uses as well. |
| `source_correction` | code | a fault of the source row this version corrects before placing the record. Codes: `nan_restored`, `accession_na`, `leucine_restored`. |
| `evidence_as_reported` | string | the evidence as the source gives it, when source_correction changed it. |
| `localised_by` | code | what localizes the methyl group to one residue, empty when nothing does (the record is then not placed). Codes: `localisation_score`, `one_methylatable_residue`, `every_arginine_dimethylated`, `curated_experimental`, `curated_database`. |
| `not_localised` | code | for a record of a study read from its own tables and not localized: why. Codes: `below_threshold`, `several_residues`, `search_disagrees`, `no_localising_evidence`, `window_fits_several`. |
| `candidate_sites` | list&lt;string&gt; | for a record that adds no site (not localized, or localized but not on one site of a protein it names): where its methyl group could be (at most 12; supporting/lower_confidence_evidence). |
| `reported_site` | string | the arginine the record reports, on the reference sequences: its accession (or an entry it names), its position and its evidence as the placement reads them, before the tests; empty when it names no arginine of an entry the record names. |
| `test_evidence_at_site` | code | site test evidence_at_site, at reported_site. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_named_entry` | code | site test named_entry. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_unique_in_proteome` | code | site test unique_in_proteome. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_methyl_on_arginine` | code | site test methyl_on_arginine. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `test_localised` | code | site test localized. Codes: `pass`, `fail`, `untested`, `curated`, `paralogs`. |
| `placed_site` | string | the site this version places the record on, empty when not placed. |
| `placement_code` | code | how the record was placed (scripts/02 remap_status). Codes: see `data_dictionary.tsv`. |
| `evidence_check` | code | independent check of the placement against the record's own evidence (scripts/qc_remap_evidence_check.py). Codes: `agrees`, `agrees_reported_arginine`, `agrees_position_only`, `agrees_no_site`, `one_of_several`, `ambiguous`, `dropped_evidence_has_site`, `contradicts`. |
| `accepted_exception` | string | why a contradiction is accepted for now (config/remap_evidence_allowlist.tsv). |
| `pending_site` | string | the site a one-sided window record supports (flag one_sided_window), empty otherwise. |
| `previous_site` | string | the site the previous version placed the record on, in its numbering. |
| `change_since_previous_version` | code | how the placement changed. Codes: `same`, `moved`, `newly_placed`, `no_longer_placed`, `not_placed`, `new_record`. |
| `reason_code` | code | why it changed. Codes: see `data_dictionary.tsv`. |

#### `supporting/lower_confidence_evidence_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `record_id` | integer | the record (supporting/source_records). |
| `study` | code | the study or database that reported it. Codes: see `data_dictionary.tsv`. |
| `reason` | code | why it adds no site. Codes: `below_threshold`, `several_residues`, `search_disagrees`, `no_localising_evidence`, `window_fits_several`, `several_sites`, `other_protein`, `not_in_reference`. |
| `candidate_sites` | list&lt;string&gt; | where the methyl group could be, on the reference sequences (at most 12). |
| `candidates_in_release` | list&lt;string&gt; | the candidates that are sites of this release (localized by another record). |
| `reported_accession` | list&lt;string&gt; | accession(s) as the source gives them. |
| `reported_gene` | string | gene symbol as the source gives it. |
| `reported_position` | string | position as the source gives it, on the accession it reports. |
| `reported_methyl_state` | string | methylation state as reported, in the codes of the site table. |
| `reported_detail` | string | what the study gives besides (localization probability, ModScore or ModLS, tables). |
| `evidence_kind` | code | the sequence evidence. Codes: `window_13`, `window_15`, `window_31`, `peptide`, `marked_peptide`, `position_only`. |
| `evidence_sequence` | list&lt;string&gt; | the peptide (marked) or window. |

#### `supporting/site_changes_v1.5.0_to_v1.6.0.tsv`

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
| `reason_codes` | list&lt;code&gt; | why it was removed, added or renumbered, or, for a kept site, why records moved onto or off it. Codes: see `data_dictionary.tsv`. |
| `records_change` | code | for a kept site, whether its records changed. Codes: `same`, `records_added`, `changed`, `replaced`. |
| `records_now_at_or_previously_at` | list&lt;string&gt; | for a removed site, where its records are placed now; for an added site, where they were placed before; for a kept site, both, for the records that moved. |
| `record_ids` | list&lt;integer&gt; | the records involved (supporting/source_records): those placed on a removed site before, or on an added site now, and for a kept site those that moved onto or off it. |

#### `supporting/published_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `version` | string | the release that had the site. |
| `published_site_id` | string | the site as that release gave it. |
| `site_ids` | list&lt;string&gt; | the same arginine in this release's numbering (several when its records now point at several), empty when none. |
| `status` | code | what it is in this release. Codes: `site`, `lower_confidence`, `not_in_reference`. |
| `tests_failed` | list&lt;code&gt; | for a lower-confidence site, the tests its best record fails. Codes: `evidence_at_site`, `named_entry`, `unique_in_proteome`, `methyl_on_arginine`, `localised`. |

#### `supporting/maron2021_detection_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `maron_table` | code | the table. Codes: `S2`, `S3`. |
| `antibody` | code | the immunoprecipitation (Table S3 only). Codes: `MMA`, `ADMA`, `SDMA`. |
| `methyl_states` | list&lt;code&gt; | the degrees called (Table S2 only). Codes: `Rme1`, `Rme2`, `Rme2a`, `Rme2s`, `Rme`. |
| `n_peptides` | integer | distinct peptides (sequence and modifications) carrying the site. |

#### `supporting/maron2021_inhibitors_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `antibody` | code | the immunoprecipitation. Codes: `MMA`, `ADMA`, `SDMA`. |
| `inhibitor` | code | the treatment, compared with DMSO. Codes: `GSK591`, `MS023`. |
| `n_peptides` | integer | peptide rows carrying the site with a measured difference. |
| `n_decreased` | integer | of these, lower with the inhibitor at adjusted p <= 0.05. |
| `n_increased` | integer | higher with the inhibitor at adjusted p <= 0.05. |
| `n_multi_site_peptides` | integer | of the peptides, those carrying another methylarginine too: their change cannot be given to one site. |
| `n_imputed_peptides` | integer | of the peptides, those with imputed intensities (Maron et al. 2021 imputed missing values; Table S3 column imputed). |
| `n_one_condition_peptides` | integer | of the peptides, those detected in one condition only (the other's three values are the imputed constant): read as detected with DMSO and not with the inhibitor (or the reverse); their p value reflects the imputation, not a measured change. |
| `log2_diff_median` | number | median log2 intensity difference over the peptides, inhibitor minus DMSO: descriptive only, not the effect of a test (it can have the sign opposite to call). |
| `log2_diff_strongest` | number | the difference of the peptide with the lowest adjusted p. |
| `adj_p_min` | number | the lowest adjusted p value (Maron et al. 2021). |
| `call` | code | over the peptides carrying the site by the shared rule (scripts/peptide_calls.py; the analyses of scripts/16 use the same): peptide-level evidence, as a peptide with several methylarginines gives each the same values. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `call_single_site_peptides` | code | the same over the peptides carrying no other methylarginine, the evidence that can be read as a change at this site; empty when every peptide carries another. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |

#### `supporting/carm1_ko_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `cell_line` | code | the cell line. Codes: `MCF7`, `MDA-MB-231`. |
| `n_peptides` | integer | peptides carrying the site. |
| `log2_ko_over_wt` | number | median log2 knockout / parental intensity over the peptides. |
| `p` | number | the lowest p value over the peptides (Shishkova et al. 2017). |
| `q` | number | the lowest BH q value over the peptides (their p values adjusted over the quantified peptides of the cell line). |
| `n_decreased` | integer | peptides at log2 <= -1 and q < 0.05. |
| `n_increased` | integer | peptides at log2 >= 1 and q < 0.05. |
| `n_multi_site_peptides` | integer | of the peptides, those carrying another methylarginine too. |
| `call` | code | a peptide decreases at log2 <= -1 and q < 0.05, increases at log2 >= 1 and q < 0.05, combined over the peptides carrying the site by the shared rule (scripts/peptide_calls.py): peptide-level evidence. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `call_single_site_peptides` | code | the same over the peptides carrying no other methylarginine; empty when every peptide carries another. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |
| `call_published` | code | the published call: twofold at nominal p < 0.05 (Shishkova et al. 2017). Codes: `decreased`, `increased`, `mixed`, `no_significant_change`. |

#### `supporting/methylarginine_site_interactions_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table). |
| `evidence_type` | code | the kind of evidence. Codes: see `data_dictionary.tsv`. |
| `localisation` | code | where the evidence falls relative to the site. Codes: `at_residue`, `candidate_residue`, `within_3`, `in_fragment`, `in_peptide`, `in_region`. |
| `partner_kind` | code | what the arginine contacts, empty when the record names no partner. Codes: `protein`, `self`, `rna`, `dna`, `other`. |
| `partner_accession` | string | the partner protein (reviewed UniProtKB accession). |
| `partner_gene` | string | its primary gene name. |
| `source` | string | PDBe-KB, ProtVar, AFDB (AlphaFold DB), IntAct, RBDmap, CAPRI, Trendel2025, UniProt, or the RNA crosslink dataset (Kramer2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS, pepR-MS_focused). |
| `structures` | list&lt;string&gt; | PDB entries (PDB rows), or the predicted model (ProtVar pair or AlphaFold DB model id). |
| `score` | number | pDockQ (ProtVar) or ipTM (AlphaFold DB) of the predicted complex. |
| `distance_a` | number | minimum interatomic distance between the arginine and the partner in the predicted model, in Å. |
| `detail` | string | for an RNA crosslink, the crosslinked peptides of that dataset with the crosslinked residue in brackets, its position and the RNA adduct masses (pepR-MS compendium Table S2); for an RBDmap fragment or a CAPRI peptide, its range and sequence; for a DNA crosslink, the crosslinked peptide with the candidate crosslinked residues in lower case (Trendel et al. 2025); bonds, pLDDT and partner residues of a predicted contact; effect, residues, PubMed and IntAct ids of an IntAct record. |
| `counts_as_evidence` | boolean | the record is site-level evidence: crosslinks at or as a candidate residue, RBDmap fragments, CAPRI peptides, PDB interfaces, predicted interfaces, confident predicted contacts, and IntAct mutations that change an interaction, short binding regions and methylated residues. interaction_status counts all of these except IntAct. |

#### `supporting/methylarginine_residue_evidence_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table. |
| `start` | integer | first residue of a run of consecutive residues with the same evidence. |
| `end` | integer | last residue. |
| `evidence_type` | code | the kind of evidence. Codes: see `data_dictionary.tsv`. |
| `localisation` | code | at the residues, or a region or peptide that holds them. Codes: `at_residue`, `candidate_residue`, `within_3`, `in_fragment`, `in_peptide`, `in_region`. |
| `n_partners` | integer | partners: most along the run (PDB protein, predicted interface) or records (IntAct). |
| `best_score` | number | best pDockQ along the run (predicted interfaces). |
| `sources` | list&lt;string&gt; | sources of the evidence. |

#### `supporting/protein_partners_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table. |
| `partner_accession` | string | STRING 12.5 physical partner at score >= 0.7 (reviewed UniProtKB accession). |
| `partner_gene` | string | its primary gene name. |
| `string_score` | number | STRING physical interaction score (0-1). |
| `partner_methylated` | boolean | the partner is a protein of the site table. |
| `partner_classes` | list&lt;code&gt; | reader or writer classes of the partner. Codes: `prmt`, `tudor`, `wd40`, `phd`. |

#### `supporting/complexes_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `complex_id` | string | Complex Portal id. |
| `complex_name` | string | its name (empty for most predicted complexes). |
| `curation` | code | curated or predicted. Codes: `curated`, `predicted`. |
| `n_members` | integer | reviewed member proteins (each protein that can fill a member slot counts). |
| `n_members_methylated` | integer | members that are proteins of the site table. |
| `methylated_members` | list&lt;string&gt; | those proteins. |
| `fold_enrichment` | number | clusterProfiler FoldEnrichment: methylated share of the members over the methylated share of reviewed human proteins. |
| `p` | number | clusterProfiler enricher p value (one-sided hypergeometric, P(X >= n_members_methylated)), against all reviewed human proteins. |
| `q` | number | Benjamini-Hochberg over the tested complexes of the same curation (complexes of at least two reviewed proteins with a methylated member). Only complexes with a methylated member are listed; those of a single reviewed protein are listed without a test. |

#### `supporting/methylarginine_nucleic_acid_contacts_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table), at a PDBe-KB RNA or DNA interface or with a nucleotide within 4 Å of its side chain in a structure analyzed. |
| `pdb_id` | string | the PDB entry (every entry PDBe-KB lists for the site with RNA or DNA is analyzed). |
| `chain` | string | the RNA or DNA chain (author chain id). |
| `nucleic_acid` | code | molecule type of the chain (RNA or DNA). Codes: `rna`, `dna`, `hybrid`. |
| `nucleic_acid_type` | code | its type, from its name, organism and preparation (DNA: nucleosomal in an entry with core histones). Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `molecule_name` | string | the chain's molecule name in the PDB entry. |
| `source_organism` | string | its source organism (synthetic construct for a made oligonucleotide). |
| `sample_preparation` | string | how the chain was made, as the entry records it: Synthetically obtained (a designed oligonucleotide), Natural source or Genetically manipulated. |
| `chain_length` | integer | nucleotides in the chain's SEQRES sequence. |
| `closest_nucleotide` | string | the nucleotide closest to the arginine side chain, numbered along the SEQRES sequence (modified nucleotides as their parent). |
| `closest_distance_a` | number | its minimum heavy-atom distance to the arginine side chain (CD, NE, CZ, NH1, NH2), in Å. |
| `contacted_nucleotides` | list&lt;string&gt; | every nucleotide of the chain within 4 Å of the arginine side chain, SEQRES numbering. |
| `contacted_author_numbers` | list&lt;string&gt; | the same nucleotides with the entry's author numbers, in the same order. |
| `contacted_distances_a` | list&lt;string&gt; | their distances in angstroms, in the same order. |
| `window_start` | integer | SEQRES number of the window's first position (below 1 when the window runs past the chain's start). |
| `window` | string | the sequence window: 8 nucleotides on each side of the nucleotide nearest to the side chain, which is the 9th; - beyond the chain's ends. Windows of different sites align on their middle. |
| `window_marked` | string | the window with the contacted nucleotides in brackets. |

#### `supporting/arginine_nucleic_acid_distances_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table with a PDB entry that PDBe-KB lists with RNA or DNA. |
| `position` | integer | an arginine of the protein modeled in at least one such entry (first model, up to two copies of the protein per entry). |
| `site_id` | string | the site, when the arginine is one (site table). |
| `methylated` | boolean | the arginine is a site of this release. |
| `entries_with_rna` | integer | entries with RNA that model the arginine. |
| `rna_distance_a` | number | the minimum heavy-atom distance from its side chain (CD, NE, CZ, NH1, NH2) to RNA over those entries, in Å, empty when never within 8 Å. |
| `rna_type` | code | the type of that RNA chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `rna_example_pdb` | string | an entry with that minimum distance. |
| `entries_with_dna` | integer | entries with DNA that model the arginine. |
| `dna_distance_a` | number | the minimum heavy-atom distance from its side chain to DNA, as for RNA. |
| `dna_type` | code | the type of that DNA chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `dna_example_pdb` | string | an entry with that minimum distance. |

#### `supporting/methylarginine_nucleic_acid_proximal_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table): a methylarginine with no side-chain atom within 4 Å of RNA or DNA at a listed interface of this complex (in none of its entries), and one within 8 Å in at least one of its entries. |
| `nucleic_acid` | code | molecule type of the nearest chain (RNA-DNA hybrids counted as RNA). Codes: `rna`, `dna`, `hybrid`. |
| `nucleic_acid_type` | code | the type of the chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `distance_a` | number | the minimum heavy-atom distance from the arginine side chain (CD, NE, CZ, NH1, NH2) to that kind of nucleic acid over the complex's entries, in Å (up to two copies of the protein per entry); 4 Å or less where PDBe-KB does not list the residue at the interface, counted as direct contact in the site table. |
| `pdb_id` | string | the entry with that minimum distance. |
| `n_entries` | integer | the entries of the complex in which a side-chain heavy atom lies within 8 Å. |
| `chain` | string | the RNA or DNA chain closest to the side chain, of this kind. |
| `molecule_name` | string | its molecule name. |
| `source_organism` | string | its source organism. |
| `sample_preparation` | string | how it was made (Synthetically obtained: a designed oligonucleotide). |
| `closest_nucleotide` | string | the nucleotide closest to the side chain, SEQRES numbering. |
| `window_start` | integer | SEQRES number of the window's first position. |
| `window` | string | 8 nucleotides on each side of the nucleotide nearest to the side chain (the 9th); - beyond the chain's ends. Aligns with the windows of supporting/nucleic_acid_contacts. |
| `window_marked` | string | the window with every nucleotide within 8 Å of the side chain in brackets. |

#### `supporting/methylarginine_nucleic_acid_near_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table): a methylarginine with no side-chain atom within 8 Å of RNA or DNA in this complex (none within 4 Å in any of its entries), and a residue at i ± 1-5 with one within 4 Å. |
| `nucleic_acid` | code | molecule type of the chain that residue contacts (RNA-DNA hybrids counted as RNA). Codes: `rna`, `dna`, `hybrid`. |
| `nucleic_acid_type` | code | the type of the chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `residue_offset` | integer | the closest such residue, in positions from the site (negative: N-terminal). |
| `residue` | integer | its position. |
| `pdb_id` | string | the entry analyzed for the complex (one per distinct set of bound sequences). |
| `n_entries` | integer | the entries with the protein and the same bound sequences. |
| `chain` | string | the RNA or DNA chain with most contacts from the residues at i ± 1-5. |
| `molecule_name` | string | its molecule name. |
| `source_organism` | string | its source organism. |
| `sample_preparation` | string | how it was made (Synthetically obtained: a designed oligonucleotide). |
| `centre_nucleotide` | string | the middle nucleotide contacted by the closest such residue, SEQRES numbering. |
| `window_start` | integer | SEQRES number of the window's first position. |
| `window` | string | 8 nucleotides on each side of the center nucleotide (the 9th); - beyond the chain's ends. Aligns with the windows of supporting/nucleic_acid_contacts. |
| `window_marked` | string | the window with every nucleotide contacted by a residue at i ± 1-5 in brackets. |

#### `supporting/nucleic_acid_regions_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein of the site table. |
| `start` | integer | first residue of a run of residues within 4 Å of RNA or DNA in a PDB structure (runs joined across gaps of up to 3 residues). |
| `end` | integer | last residue of the run. |
| `nucleic_acid` | code | molecule type of the chains the residues contact (RNA-DNA hybrids counted as RNA). Codes: `rna`, `dna`, `hybrid`. |
| `n_entries` | integer | PDB entries in which the run contacts it (every entry PDBe-KB lists for the protein with RNA or DNA). |
| `n_complexes` | integer | distinct complexes analyzed: one entry per set of bound sequences, the others counted with it. |
| `methylarginines` | list&lt;string&gt; | the sites (site table) inside the run. |
| `bound_segments` | list&lt;string&gt; | the stretches of RNA or DNA the run contacts, 5' to 3': contacted nucleotides in brackets with 2 more each side, and … between clusters further apart than 4 nucleotides; at most 6, most entries first. |
| `segment_entries` | list&lt;string&gt; | PDB entries showing each segment, in the same order. |
| `segment_example_entries` | list&lt;string&gt; | an entry showing each segment. |
| `segment_molecule_names` | list&lt;string&gt; | the molecule name of each segment's chain. |
| `segment_types` | list&lt;code&gt; | the type of each segment's chain. Codes: `rrna`, `trna`, `snrna`, `mrna`, `other_ncrna`, `viral_rna`, `rna_oligo`, `other_rna`, `nucleosomal_dna`, `telomeric_dna`, `other_dna`. |
| `n_segments` | integer | distinct segments, including those not listed. |

#### `supporting/rna_bind_n_seq_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | methylated protein with an ENCODE RNA Bind-n-Seq experiment (RNA sequences enriched by the purified protein from a random pool in vitro). |
| `encode_experiment` | string | the ENCODE experiment. |
| `concentration` | string | the protein concentration at which its top 5-mer is most enriched; the 5-mers are ranked there. |
| `top_5mers` | list&lt;string&gt; | the 12 most enriched 5-mers. |
| `top_5mer_R` | list&lt;string&gt; | their enrichment R (frequency in the bound RNA over the input), in the same order. |
| `logo_consensus` | list&lt;string&gt; | the consensus (IUPAC) of each logo ENCODE gives for the experiment, strongest first. |

#### `supporting/dna_motifs_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | methylated protein with a JASPAR 2026 CORE vertebrates non-redundant profile (its own, or a heterodimer's it is part of). |
| `matrix_id` | string | the JASPAR matrix (jaspar.elixir.no/matrix/<id>). |
| `name` | string | the matrix name (A::B for a heterodimer). |
| `tf_class` | string | the transcription factor class, as JASPAR gives it. |
| `tf_family` | string | the transcription factor family. |
| `length` | integer | positions in the matrix. |
| `consensus` | string | IUPAC consensus of the matrix. |

#### `supporting/htr_selex_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | methylated protein with an HTR-SELEX model (Jolma et al., Genome Res 30:962, 2020: RNA selected by the purified protein from random 40-nucleotide RNA in vitro, CC BY-NC 4.0). |
| `model_id` | string | the model: jolma2020_linear_<n> from Supplemental Table S2, jolma2020_structured_<n> from Table S3, in their order. |
| `motif_kind` | code | whether the motif is linear or lies within a stem-loop. Codes: `linear`, `structured`. |
| `motif_type` | string | Primary or Secondary, as the paper gives it. |
| `experiment` | string | the experiment that gave the model (Experiment_1, or its Replicate). |
| `cycle` | string | the selection cycle the model was derived from. |
| `oligomer` | string | monomeric or dimeric binding, as the paper gives it (linear models). |
| `length` | integer | positions in the model. |
| `consensus` | string | IUPAC consensus of the position frequency matrix. |

#### `supporting/eclip_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | methylated protein with an eCLIP experiment (RNA crosslinked to it in K562 or HepG2 cells). |
| `source` | string | ENCODE (reproducible peaks) or ENCORE (Her et al. 2026, GEO GSE315407: reproducible enriched windows of at most 100 nt, Skipper). |
| `experiment` | string | the ENCODE experiment (ENCSR id) or the ENCORE experiment (its GEO file name: RBP, cell line, Yeo-lab number). |
| `cell_line` | string | K562, HepG2, or another biosample. |
| `n_peaks` | integer | reproducible peaks (ENCODE, replicates 1 and 2, GRCh38) or reproducible enriched windows (ENCORE); the peaks_ columns count the same. |
| `peaks_CDS` | integer | peaks whose middle falls in a coding sequence of a GENCODE v50 basic transcript of their strand. |
| `peaks_UTR3` | integer | in a 3' UTR (and in no coding sequence). |
| `peaks_UTR5` | integer | in a 5' UTR. |
| `peaks_noncoding_exon` | integer | in an exon of a non-coding transcript. |
| `peaks_intron` | integer | in an intron. |
| `peaks_intergenic` | integer | outside every basic transcript. |
| `top_target_genes` | list&lt;string&gt; | the 25 genes with most peaks. |
| `top_target_peaks` | list&lt;string&gt; | their peaks (or windows), in the same order. |
| `top_5mers` | list&lt;string&gt; | the 5 5-mers (RNA alphabet) most enriched in the top 1,000 peaks or windows against those of the other RBPs of the same source, weighted to the same mix of regions (scripts/43e; ranked by a two-proportion z, a ranking only). |
| `top_5mer_ratios` | list&lt;string&gt; | their frequency ratios, in the same order. |
| `same_data_as` | string | ENCORE rows: the released ENCODE experiment whose data this re-analyses, if any (not independent evidence). |

#### `supporting/disease_terms_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `term_id` | string | DisGeNET term (UMLS CUI). |
| `term_name` | string | its name. |
| `family` | string | its disease family (scripts/disease_families.py): audited, or else additional; empty when none. |
| `n_genes_term` | integer | genes of the term in the universe. |
| `n_methylated` | integer | of these, the genes of methylarginine-containing proteins. |
| `fold_enrichment` | number | observed over expected. |
| `p` | number | hypergeometric p value (clusterProfiler enricher on the DisGeNET terms of DOSE 4.4.0). |
| `p_adjust` | number | Benjamini-Hochberg adjusted p value over the tested terms. |
| `q` | number | q value (qvalue, lambda 0.05, as DOSE). |
| `genes` | list&lt;string&gt; | the genes of the methylarginine-containing proteins in the term (symbols, org.Hs.eg.db). |
| `uniprot_accessions` | list&lt;string&gt; | the proteins of the site table encoded by those genes. |

#### `supporting/disease_terms_tested_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `term_id` | string | DisGeNET term (UMLS CUI). |
| `term_name` | string | its name (a ";" inside a name, the list separator, is written ","). |
| `family` | string | its disease family (scripts/disease_families.py); empty when none. |
| `n_genes_term` | integer | genes of the term in DOSE's background (reviewed human genes with a DisGeNET term). |
| `n_methylated` | integer | of these, the genes of methylarginine-containing proteins. |
| `fold_enrichment` | number | fold enrichment as DOSE reports it. |
| `p` | number | one-sided hypergeometric p value (clusterProfiler enricher on the DisGeNET terms of DOSE 4.4.0). |
| `p_adjust` | number | Benjamini-Hochberg adjusted p value over the tested terms. |
| `q` | number | q value (DOSE). |
| `enriched` | boolean | p, p_adjust and q at most 0.05, the cut of the published analysis (supporting/disease_terms). |
| `genes` | list&lt;string&gt; | the genes of the methylarginine-containing proteins in the term (symbols). |
| `uniprot_accessions` | list&lt;string&gt; | the proteins of the site table encoded by those genes. |

#### `supporting/disease_ontology_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `term_id` | string | Human Disease Ontology term (DOID). |
| `term_name` | string | its name (a ";" inside a name is written ","). |
| `parents` | list&lt;string&gt; | its parent terms (is_a). |
| `n_genes_term` | integer | genes of the term and of the terms under it, among the reviewed human genes with a Disease Ontology term (the gene sets of DOSE enrichDO). |
| `n_methylated` | integer | of these, the genes of methylarginine-containing proteins. |
| `fold_enrichment` | number | fold enrichment (clusterProfiler). |
| `p` | number | one-sided hypergeometric p value (clusterProfiler::enricher on the Disease Ontology gene sets; the same as DOSE enrichDO). |
| `q` | number | Benjamini-Hochberg adjusted p value over the tested terms (at least 10 genes). |
| `genes` | list&lt;string&gt; | the genes of the methylarginine-containing proteins in the term (symbols). |
| `uniprot_accessions` | list&lt;string&gt; | the proteins of the site table encoded by those genes. |

#### `supporting/disease_families_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `family` | string | disease family (scripts/disease_families.py): the eight audited families, then the additional ones, each in order of significance. |
| `audited` | boolean | one of the eight audited families (scripts/disease_families.py FAMILY_RULES, the original pipeline's hand-checked rules); otherwise an additional, rule-based family. |
| `n_terms` | integer | DisGeNET terms assigned to it. |
| `n_genes` | integer | genes of the universe in at least one of its terms. |
| `n_methylated_genes` | integer | of these, the genes of methylarginine-containing proteins. |
| `n_genes_tested` | integer | the background of the test: reviewed human genes with at least one DisGeNET term (the same for every family). |
| `n_methylated_genes_tested` | integer | of these, the genes of methylarginine-containing proteins (the same for every family). |
| `fold_enrichment` | number | fold enrichment of the family (clusterProfiler::enricher against the reviewed human genes with a DisGeNET term; statistics: disease_family). |
| `q` | number | its Benjamini-Hochberg adjusted p value (over the audited families, or over the additional ones). |

#### `supporting/arginine_background_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | a protein whose sites the residue test counts. |
| `analysis` | string | the residue test (statistics analysis): motif, disorder, domain_context, domain_class or structure_context. |
| `category` | string | the class or motif of the test; * for all of the protein's other arginines. |
| `other_arginines` | integer | the arginines of the protein that are not sites of the release, in the class or with the motif: the background of the test, which compares the sites with the other arginines of the same proteins. |

#### `supporting/density_windows_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `window_aa` | integer | window width in residues. |
| `start` | integer | first residue of the densest window of that width. |
| `end` | integer | last residue. |
| `sites` | integer | methylarginine sites in it. |
| `arginines` | integer | arginines in it. |
| `null_mean_sites` | number | mean, under the null, of the maximum number of sites in any window of that width (sites placed on the protein's arginines at random). |
| `log2_obs_over_null` | number | log2 of observed over null mean. |
| `p` | number | empirical p value (1,000 permutations). |
| `q` | number | Benjamini-Hochberg q value over all protein-window tests. Only windows with q <= 0.05 are listed. |
| `site_positions` | list&lt;integer&gt; | the sites in the window. |

#### `supporting/condensates_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `condensate_id` | string | CD-CODE condensate id. |
| `condensate_name` | string | its name. |
| `condensate_type` | string | its type in CD-CODE. |
| `n_members` | integer | human proteins CD-CODE lists as members. |
| `n_members_methylated` | integer | of these, proteins of the site table. |

#### `supporting/hpa_locations_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `gene` | string | gene as HPA gives it. |
| `ensembl_gene` | string | Ensembl gene HPA maps to the protein. |
| `location` | string | subcellular location. |
| `main` | boolean | a main location (otherwise additional). |
| `reliability` | string | HPA reliability: Enhanced, Supported, Approved or Uncertain. |

#### `supporting/other_ptm_sites_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | protein. |
| `position` | integer | position on the reference sequence. |
| `residue` | string | the modified residue. |
| `ptm_type` | code | the modification. Codes: `phospho`, `lys_acyl`, `ub_sumo`, `lys_methyl`. |
| `source` | string | ptm_master: the iPTMnet-format PTM master table (rows on the canonical accession, residue checked against the reference sequence). |
| `evidence` | list&lt;string&gt; | the databases the rows come from, as the master table names them. |
| `n_records` | integer | rows of the master table. |

#### `supporting/writer_perturbations_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site. |
| `perturbation` | code | the knockout, knockdown or inhibitor, and the cells. Codes: `carm1_ko_mcf7`, `carm1_ko_mdamb231`, `prmt7_ko_hct116`, `prmt1_kd_hek293t`, `prmt5i_epz015666_hela`, `prmt7_kd_hek293`, `prmt7i_sgc3027_hek293`, `carm1_kd_hek293`, `prmt5_kd_hek293`. |
| `n_peptides` | integer | quantified peptides carrying the site. |
| `log2_treated_over_control` | number | median log2 perturbed / control over those peptides: descriptive, not the effect of a test. |
| `p` | number | the lowest p value over the peptides, as the study gives it. |
| `q` | number | the lowest q over the peptides: the study's adjusted p where it gives one, else its p adjusted here (Benjamini-Hochberg, within the comparison); empty when the study gives no p. |
| `n_decreased` | integer | peptides called decreased by the criteria of `call`. |
| `n_increased` | integer | peptides called increased by the criteria of `call`. |
| `n_multi_site_peptides` | integer | of the peptides, those carrying another methylarginine too. |
| `call` | code | the call: `call_twofold_q05` or `call_published`, as `call_criteria` says (the perturbation's code says why); peptide-level evidence. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `call_criteria` | code | which of the two calls `call` is. Codes: `twofold_q05`, `published`. |
| `call_single_site_peptides` | code | the same over the peptides that carry no other methylarginine and whose one mark the search localized; empty when there are none. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `call_twofold_q05` | code | the same criteria for every perturbation: a peptide changes at |log2| >= 1 and q < 0.05, combined over the peptides carrying the site (scripts/peptide_calls.py); not_tested where the study gives no p value. Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |
| `call_published` | code | by the published criteria (code perturbation); empty where the study makes no call (as for a dimethylated site in a monomethylation list). Codes: `decreased`, `increased`, `mixed`, `no_significant_change`, `not_tested`. |

#### `supporting/functional_experiments_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `experiment_id` | string | FE and four digits, never reused. |
| `site_id` | string | the site, accession:R position; empty when the experiment is not resolved to one site. |
| `in_release` | boolean | the site is a site of this release. |
| `uniprot_accession` | string | reviewed human accession. |
| `gene` | string | its primary gene name. |
| `position` | integer | the arginine; empty when the experiment is not resolved to one site. |
| `region` | string | for a row without a site, which arginines the experiment concerns. |
| `species` | string | of the cells or protein in the experiment. |
| `system` | string | cells or organism, or the in vitro system. |
| `writer` | string | the PRMT the paper names as methylating the site, empty when none. |
| `experiment_kind` | code | what kind of experiment; empty when only the abstract was read. Codes: `methylated_vs_unmethylated`, `site_mutation`, `prmt_inhibitor`, `prmt_knockdown`, `prmt_knockout`, `correlation`, `prediction`. |
| `perturbation` | string | what was done, as the paper names it. |
| `perturbation_note` | string | what else the perturbation changes. |
| `readout` | string | what was measured. |
| `direction` | code | of the readout under the perturbation. Codes: `increased`, `decreased`, `no_change`. |
| `finding` | string | one sentence, close to the paper's words. |
| `figure` | string | where in the paper. |
| `reference` | string | the paper. |
| `doi` | string | its DOI. |
| `pmid` | string | its PubMed id. |
| `read_from` | code | what the curator read. Codes: `full_text`, `abstract`. |
| `curated_on` | string | date of curation. |
| `note` | string | anything else. |

#### `supporting/eclip_rna_compartments_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `uniprot_accession` | string | reviewed accession of the protein with eCLIP data. |
| `protein` | string | its name in the eCLIP experiments. |
| `methylated` | boolean | the protein has a site in this release. |
| `map` | code | the map of RNA or DNA. Codes: `apexseq`, `nuclear_body`, `tsaseq`, `hic`. |
| `compartment` | string | the compartment (a subcompartment for Hi-C). |
| `bound_rnas_placed` | integer | the protein's bound RNAs (genes) the map places. |
| `bound_rnas_in_compartment` | integer | of these, those the map places in the compartment. |
| `odds_ratio` | number | Mantel-Haenszel odds ratio of its bound RNAs being in the compartment against the other RNAs the map places, within deciles of expression in the map's cells crossed with deciles of gene length (GENCODE v19; the genes v19 lacks are one more length stratum within each expression decile). |
| `odds_ratio_ci_low` | number | lower bound of its 95% confidence interval (Robins-Breslow-Greenland). |
| `odds_ratio_ci_high` | number | upper bound. |
| `p` | number | continuity-corrected Cochran-Mantel-Haenszel test, two-sided. |
| `q` | number | Benjamini-Hochberg over every protein and compartment of the map. |
| `gathers` | boolean | its bound RNAs gather in the compartment: q <= 0.05 and odds ratio above 1. |

#### `supporting/methylarginine_site_interaction_records_v1.6.0.tsv`

| Column | Type | Description |
|---|---|---|
| `site_id` | string | the site (site table). |
| `record_class` | code | what the record shows the partner doing at the site. Codes: see `data_dictionary.tsv`. |
| `evidence_types` | list&lt;code&gt; | the rows of supporting/site_interactions it comes from (pdb_methylated_contact: an entry modeling the site methylated that PDBe-KB does not list for it), in code order. Codes: see `data_dictionary.tsv`. |
| `partner_kind` | code | what the partner is. Codes: `protein`, `self`, `rna`, `dna`, `peptide`, `gene`, `small_molecule`, `other`. |
| `partner_accession` | string | its UniProtKB accession, of any species, empty when it has none. |
| `partner_gene` | string | its gene (human reviewed entries), or the name IntAct gives. |
| `partner_name` | string | its name in the experiment: the PDB entity description or the IntAct participant description. |
| `partner_species` | string | the partner's organism in the experiment. |
| `source` | string | PDB, IntAct, or the study or model (CAPRI, pCLAP, RBDmap, RBR-ID, Trendel2025, the pepR-MS datasets, UniProt, ProtVar, AFDB). |
| `source_id` | string | the PDB entry, the IntAct interaction (EBI-), or the predicted model. |
| `method` | string | how it was found: the PDB experimental method and resolution, the IntAct detection method (PSI-MI), or the study's method. |
| `context` | string | the PDB entry title; for IntAct the interaction type and host organism; the localization of a crosslink. |
| `site_state` | code | the site in that experiment. Codes: `arginine`, `Rme1`, `Rme2a`, `Rme2s`, `Rme2`, `Rme`, `other_residue`, `mutated`, `not_read`, `not_recorded`. |
| `state_detail` | string | the residue modeled (chemical component, how the entity was matched, the copies), the PSI-MOD term or the mutation. |
| `other_marks` | list&lt;string&gt; | other modified residues of the site's protein in the experiment (K9me3, S10ph; UniProt numbering where the entity aligns). |
| `site_species` | string | the organism of the site's protein in the experiment. |
| `detail` | string | the partner residues within 4 Å of a methylated site, the IntAct features of both participants, or the evidence detail of supporting/site_interactions. |
| `pubmed` | string | PubMed id of the paper. |
| `doi` | string | its DOI. |
| `citation` | string | first author, journal and year (the primary citation of a PDB entry). |
| `paper_title` | string | the title of the paper. |
| `flags` | list&lt;code&gt; | what to check, in code order. Codes: `partner_other_species`, `site_protein_other_species`, `matched_by_sequence`, `state_not_reported`, `paper_names_other_arginine`, `title_names_paralogue`, `state_unexpected_for_enzyme`, `no_publication`, `negative`. |

## Checksums

Run `md5sum -c MD5SUMS` in this folder, and in `analyses/` and `analyses_with_lower_confidence/`.
