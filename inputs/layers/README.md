# Layer inputs

These tables hold the annotation layers. The website shows them on protein and site pages; none of them adds
methylarginine sites to the union. Each layer is registered in `scripts/release_layers.py`, which lists its files (the
release build checks them against `../MANIFEST.tsv` first) and the columns, tables and tests each adds to
the release (`../../ADDING_A_SOURCE.md` says how to add one).

| Layer | Made by | Files |
|---|---|---|
| CPLM 4.0 lysine methylation | `scripts/39_stage_cplm4_lys_methylation.py` | `cplm4_human_lys_methyl_validated.tsv`, `cplm4_human_lys_methyl_rejected.tsv`, `cplm4_summary.json` |
| CARM1 knockout response (Shishkova et al. 2017) | `scripts/41_stage_carm1_ko_sites.py` | `carm1_ko_peptides.tsv`, `carm1_ko_sites.tsv`, `carm1_ko_summary.json` |
| Human Protein Atlas 25.1 subcellular locations | `scripts/42_stage_hpa_subcellular.py` | `hpa_subcellular_locations.tsv`, `hpa_summary.json` |
| Cancer cell lines and tumours: DepMap Public 24Q4 (CRISPR gene effect, changes at arginines, MTAP copy number), PRISM Repurposing 24Q2 (PRMT inhibitors) and the TCGA MC3 somatic calls at arginines (GDC open access) | `scripts/44b_extract_tcga_mc3.py` (from the MC3 public MAF, the GDC cases of TCGA and UniProt's Ensembl cross-references), then `scripts/44_stage_depmap.py`; the release build adds its columns and tests (`scripts/release_layers.py`) and `scripts/44c_export_depmap.py` writes the website's `docs/data/depmap/` from the release | `depmap_24q4_gene_effect_summary.tsv`, `depmap_24q4_arginine_changes.tsv`, `depmap_24q4_models.tsv`, `depmap_24q4_prmt_gene_effect.tsv`, `prism_24q2_prmt_inhibitors.tsv`, `tcga_mc3_arginine_changes.tsv`, `tcga_mc3_projects.tsv`, `depmap_summary.json` |
| Cancer changes at every candidate K, S, T and Y (v1.6.0), and each protein's background: the TCGA MC3 calls and DepMap 24Q4 by the rules of the arginine layer | `scripts/66_stage_cancer_changes.py` (from the MC3 public MAF and OmicsSomaticMutations.csv of the arginine layer, the GDC cases of TCGA, UniProt's Ensembl cross-references and the genes DepMap screened) | `tcga_mc3_changes_{K,S,T,Y}.tsv.gz`, `depmap_24q4_changes_{K,S,T,Y}.tsv.gz`, `cancer_change_protein_background.tsv.gz`, `cancer_changes_summary.json` |
| ClinVar variants at every candidate R, K, S, T and Y (v1.6.0), and each protein's background | `scripts/67_stage_clinvar.py` (from ClinVar's hgvs4variation and variant_summary tables and UniProt's RefSeq cross-references) | `clinvar_changes_{R,K,S,T,Y}.tsv.gz`, `clinvar_protein_background.tsv.gz`, `clinvar_summary.json` |
| gnomAD v4.1 exome variants at every candidate R, K, S, T and Y (v1.6.0), and each protein's background | `scripts/69_stage_gnomad.py` (streamed from Google's public copy, `gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/`) | `gnomad_v4_1_changes_{R,K,S,T,Y}.tsv.gz`, `gnomad_protein_background.tsv.gz`, `gnomad_summary.json` |
| Writer perturbations: PRMT7 knockout (Szewczyk et al. 2020), PRMT1 knockdown (Hartel et al. 2019), PRMT5 inhibitor (Lu et al. 2024), PRMT7, CARM1 and PRMT5 knockdowns and the PRMT7 inhibitor SGC3027 (Li et al. 2021) | `scripts/45_stage_prmt7_szewczyk2020.py`, `scripts/46_stage_hartel2019_prmt1_kd.py`, `scripts/51_stage_lu2024_prmt5.py`, `scripts/50_stage_li2021_perturbations.py` | `prmt7_ko_*`, `hartel2019_prmt1_kd_*`, `lu2024_prmt5_*`, `li2021_*` |
| PaxDb 6.1 protein abundance | `scripts/47_stage_paxdb.py` | `paxdb_human_integrated_abundance.tsv`, `paxdb_summary.json` |
| Reactome 97 pathways | `scripts/48_stage_reactome.py` | `reactome_human_protein_pathways.tsv`, `reactome_human_pathways.tsv`, `reactome_summary.json` |

- `../MANIFEST.tsv` (one manifest for all of `inputs/`) describes each file: rows, MD5, what it holds, source, version, citation, licence and how it
  was made.
- `raw_files.tsv` lists the downloads the tables were made from, with URL, size, MD5, SHA-256 and download date.
- The summary JSON files are the scripts' own summaries, with paths shortened to file names. Their input
  checksums equal those stored in the website data (`docs/data/methylarg_sites.json`, `meta.cplm4`, `meta.hpa`,
  `meta.carm1_ko`).

## Checks

On 20 October 2026 the three scripts were run again from the downloads in `raw_files.tsv`. All five tables came
out byte-identical, both times:

- with the FASTA the scripts originally used, the UniProt 2026_03 reviewed human reference proteome
  (20,416 sequences). That file is now in Git as `../reference/uniprot_human_reviewed_2026_03.fasta`;
- with `../interactions/uniprot_human_reviewed_2026-10-16.fasta`, also in Git.

On 21 October 2026 the HPA tables were made again after a fix to scripts/42. Where HPA maps one accession to
several genes, the script had kept the reliability and main/additional flag of the first gene listed; it now keeps
the best reliability and counts a location as main when any of the genes has it as main. Four rows changed.

## What each table holds

- **`cplm4_human_lys_methyl_validated.tsv`** (5,356 sites on 2,939 proteins): CPLM 4.0 human lysine methylation
  sites whose CPLM sequence matches the canonical UniProt sequence. The match is either an identical sequence, or
  the same ±7-residue window around the lysine. Nothing is remapped. The CPLM id, evidence label (`Exp.`, `Dat.`,
  `Exp./Dat.`) and references are given as CPLM gives them.
- **`cplm4_human_lys_methyl_rejected.tsv`** (559 rows): the rows that failed, with the reason:
  - 494: the accession is not in the reference proteome;
  - 60: no lysine at that position in the canonical sequence;
  - 5: the window differs.
- **`carm1_ko_peptides.tsv`** (515 peptides): ADMA peptides quantified in parental and CARM1 knockout MCF7 and
  MDA-MB-231 cells. Each peptide is placed on the canonical sequence and carries its log2 KO/WT, p value and call:
  - down: log2 ≤ −1 and p < 0.05;
  - up: log2 ≥ 1 and p < 0.05.
- **`carm1_ko_sites.tsv`** (603 rows, 431 sites): one row per methylarginine and cell line, with:
  - median log2 KO/WT;
  - lowest p;
  - peptide counts;
  - the call (down, up, unchanged, mixed).
- **`hpa_subcellular_locations.tsv`** (27,150 rows on 13,301 proteins): one row per protein and HPA location,
  with:
  - whether it is a main (1) or additional (0) location;
  - the HPA reliability (Enhanced, Supported, Approved, Uncertain).

- **`tcga_mc3_arginine_changes.tsv`** (274,408 rows, 9,064 tumour samples): somatic missense and stop changes at
  arginines of reviewed canonical proteins in TCGA tumours, from the MC3 public MAF of the GDC
  (`mc3.v0.2.8.PUBLIC.maf.gz`, Ellrott et al. 2018; the calls of the PanCancer Atlas). Calls are kept as MC3 kept them:
  FILTER `PASS`, and in OV and LAML also `wga` alone. The change is read on a transcript that UniProt 2026_03 assigns to
  the entry's displayed isoform, from the MAF's `all_effects` (`mapped_by` `canonical_transcript`, 95% of rows); where
  the variant falls in no such transcript, on the MAF's own transcript, placed through the UniProt entry it names or
  the gene symbol (`maf_transcript`). Either way a change counts only where the canonical sequence has an arginine at
  its position. GDC open-access data carry no licence; the GDC's data access policies (NIH Genomic Data Sharing policy)
  ask that TCGA and the MC3 publication be acknowledged and that no one try to identify participants.
- **`tcga_mc3_projects.tsv`** (33 rows): each TCGA project, with its tumour samples in the MAF (10,294 in all), those
  with a kept call, and the arginine changes extracted.

- **`tcga_mc3_changes_{K,S,T,Y}.tsv.gz`** and **`depmap_24q4_changes_{K,S,T,Y}.tsv.gz`**: the somatic missense and stop
  changes at every lysine, serine, threonine and tyrosine of reviewed canonical proteins, a known site or not, in TCGA
  tumors and DepMap cell lines, for the lysine, methyllysine and phosphorylation site sets. They have the columns of the
  arginine tables, and `gene`, and follow their rules: rerun on arginines, `scripts/66_stage_cancer_changes.py` writes
  `tcga_mc3_arginine_changes.tsv` and `depmap_24q4_arginine_changes.tsv` byte for byte. The Ensembl transcripts of the
  canonical isoforms come from UniProt 2026_04, the current release on the day of staging; positions are checked on the
  2026_03 reference, as for arginines. TCGA rows: K 61,510, S 142,576, T 59,012, Y 25,757; DepMap rows: K 20,835,
  S 44,341, T 24,605, Y 13,533. The files are gzip-compressed without a timestamp, so restaging gives the same bytes.
  `tcga_mc3_changes_R.tsv.gz` and `depmap_24q4_changes_R.tsv.gz` are the arginine tables' rows, byte for byte, with `gene`;
  the release reads them, and the arginine tables stay as they are.
- **`gene`**, in every change table here (TCGA, DepMap, ClinVar, gnomAD): the HGNC approved symbol of the gene the
  change lies in. It matters where several genes encode one protein: the 14 genes of histone H4 (P62805) are pooled
  under one accession. The sources name genes in different years' terms, so each row's gene is resolved through a stable
  identifier where the source gives one (the Ensembl 105 gene of the transcript for TCGA and gnomAD, DepMap's Ensembl
  gene, ClinVar's NCBI Gene ID), else through its symbol (approved, a unique previous symbol, a unique alias), with the
  HGNC complete set of 13 November 2026 (CC0, https://genenames.org/about/license) and the Ensembl 105 GTF
  (`scripts/ensembl_genes.py`).
- **`cancer_change_protein_background.tsv.gz`** (122,496 rows): for each reviewed canonical entry and residue class (R, K,
  S, T, Y, and all residues), the residues of that class in its sequence and the TCGA and DepMap changes on them, with
  the tumor samples or cell lines that carry one and the positions changed. A site's changes can be compared with those
  of the other residues of its kind in the same protein. The denominators (10,294 tumor samples, 1,929 cell lines with
  mutation calls) are in `cancer_changes_summary.json`.

- **`clinvar_changes_{R,K,S,T,Y}.tsv.gz`**: the ClinVar variants (files dated 15 and 16 November 2026) that change a
  residue of the site sets on a reviewed canonical protein: a substitution or a stop written on a RefSeq protein that
  UniProt assigns to the entry's displayed isoform, kept where the 2026_03 canonical sequence has that residue. One row
  per entry, position, new residue and variation, with the germline classification, review status and ClinVar's
  stars, submitters, origin, conditions, the somatic clinical impact and oncogenicity, and the dbSNP ID; the variation's
  page (`https://www.ncbi.nlm.nih.gov/clinvar/variation/<variation_id>/`) lists its records. Rows: R 335,933, K 99,522,
  S 186,004, T 142,865, Y 64,910. NCBI places no restrictions on the use or distribution of ClinVar data; submitters
  may claim rights in what they submitted.
- **`clinvar_protein_background.tsv.gz`** (122,496 rows): for each entry and residue class, the residues of that class,
  the ClinVar variants on them, those classified pathogenic or likely pathogenic, and the positions changed.

- **`gnomad_v4_1_changes_{R,K,S,T,Y}.tsv.gz`**: the gnomAD v4.1 exome variants (730,947 exomes) that pass every gnomAD
  filter and that VEP calls missense or stop gained on a transcript UniProt assigns to a reviewed entry's displayed
  isoform, kept where the 2026_03 canonical sequence has that residue. One row per entry, position, new residue and
  variant (GRCh38 chromosome-position-reference-alternate), with the dbSNP ID, allele count and number, frequency,
  homozygotes, the ancestry group with the highest frequency and its frequency, and the highest filtering allele
  frequency (95%). Rows: R 1,133,321, K 517,573, S 965,320, T 678,065, Y 307,801. The 24 sites VCFs (198 GB) were streamed,
  not stored; `raw_files.tsv` gives the bytes and MD5 read for each. The primary data are CC0; SpliceAI and CADD, which
  carry non-commercial terms, are not read.
- **`gnomad_protein_background.tsv.gz`** (122,496 rows): for each entry and residue class, the residues of that class,
  the gnomAD variants on them, the positions changed and the summed allele count.

## CPLM tables

CPLM 4.0 states "Copyright © 2010-2026 The CUCKOO Workgroup. All Rights Reserved" and gives no licence
(cplm.biocuckoo.cn, checked 20 October 2026). For that reason (the context marks, `../marks/`, do not include CPLM since 17 November 2026):

- the two CPLM tables of this folder stay out of Git (`.gitignore`), and are marked `checksum_only` in `../MANIFEST.tsv`, which
  gives the MD5 of both;
- the CPLM summary is in Git;
- anyone who downloads `Methylation.zip` (MD5 in `raw_files.tsv`) can rebuild them with scripts/39.

The offline archive keeps them in a separate, restricted file.
