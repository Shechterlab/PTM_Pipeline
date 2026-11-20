# Changes from v1.5.1 to v1.6.0

Everything of v1.5.1 is in v1.6.0 under this release's names, its data unchanged. v1.6.0 adds the somatic changes at every candidate residue of the site sets, in TCGA tumors and DepMap cell lines, and the germline variants of ClinVar and gnomAD there.

## New: somatic changes at every candidate residue

`supporting/tcga_mc3_changes_<residue>_v1.6.0.tsv.gz` and `supporting/depmap_24q4_changes_<residue>_v1.6.0.tsv.gz` list the somatic missense and stop changes at every arginine, lysine, serine, threonine and tyrosine of the reviewed canonical proteins, a known site or not: in TCGA tumors (the MC3 calls of the PanCancer Atlas, 10,294 tumor samples in 33 projects) and in DepMap 24Q4 cell lines (1,929 lines with mutation calls). One row per tumor sample or cell line and change. The column `site_sets` names the site sets of this release with a site at that residue (empty: a candidate residue, not a known site), so a set's changes can be read at every candidate residue or at its known sites only. The column `gene` names the gene each change lies in (its current HGNC approved symbol, for all four sources), which matters where several genes encode one protein: the 14 genes of histone H4 (P62805), for example, are pooled under one accession.

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

`supporting/clinvar_changes_<residue>_v1.6.0.tsv.gz` and `supporting/gnomad_v4_1_changes_<residue>_v1.6.0.tsv.gz` list the human germline variants that change an arginine, lysine, serine, threonine or tyrosine of a reviewed canonical protein (a substitution or a stop), a known site or not, with `site_sets` as in the cancer tables. ClinVar (files of 2026-11-15 (hgvs4variation), 2026-11-16 (variant_summary)): one row per variation, with its germline classification, review status and stars, submitters, origin, conditions, somatic clinical impact and oncogenicity; the variant counts where one of its protein changes is written on a RefSeq protein UniProt assigns to the entry's displayed isoform. gnomAD v4.1 (730,947 exomes): one row per variant that passes every gnomAD filter, with its allele count and number, frequency, homozygotes, the ancestry group with the highest frequency and the highest filtering allele frequency; the change is VEP's, on a transcript UniProt assigns to the displayed isoform. Either way a variant counts only where the 2026_03 canonical sequence has that residue. As in the cancer tables, `gene` names the gene of each variant.

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

## Changes not shown

`supporting/variant_changes_excluded_v1.6.0.tsv.gz` lists the changes a source places on an entry that the change tables leave out because the residue the source names there differs from the canonical sequence's (UniProt 2026_03), or the position lies beyond it, where either residue is R, K, S, T or Y: for example a gene of a protein family whose own protein differs from the canonical one at that position. 42,056 changes on 1,240 proteins: TCGA MC3 32,295, DepMap 24Q4 611, ClinVar 2,072, gnomAD v4.1 7,078.

## Corrected: RBDmap and pCLAP in the interaction notes

The interaction notes of v1.5.1's README said, under what the tables of methyllysine, lysine and phospho lack, that the extents of RBDmap fragments and pCLAP peptides are not republished. They are: the methylarginine interaction tables carry each RBDmap fragment and pCLAP peptide with its extent and sequence, credited (RBDmap: Castello et al. 2016, CC BY 4.0; pCLAP: Mullari et al. 2017, called from PRIDE PXD005739 with the rule of Mullari et al. 2023). This README says so; v1.5.1's published folder is left as archived.

## New files

- `supporting/cancer_change_protein_background_v1.6.0.tsv.gz`
- `supporting/clinvar_changes_K_v1.6.0.tsv.gz`
- `supporting/clinvar_changes_R_v1.6.0.tsv.gz`
- `supporting/clinvar_changes_S_v1.6.0.tsv.gz`
- `supporting/clinvar_changes_T_v1.6.0.tsv.gz`
- `supporting/clinvar_changes_Y_v1.6.0.tsv.gz`
- `supporting/clinvar_protein_background_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_change_models_v1.6.0.tsv`
- `supporting/depmap_24q4_changes_K_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_R_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_S_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_T_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_Y_v1.6.0.tsv.gz`
- `supporting/gnomad_protein_background_v1.6.0.tsv.gz`
- `supporting/gnomad_v4_1_changes_K_v1.6.0.tsv.gz`
- `supporting/gnomad_v4_1_changes_R_v1.6.0.tsv.gz`
- `supporting/gnomad_v4_1_changes_S_v1.6.0.tsv.gz`
- `supporting/gnomad_v4_1_changes_T_v1.6.0.tsv.gz`
- `supporting/gnomad_v4_1_changes_Y_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_K_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_R_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_S_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_T_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_Y_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_projects_v1.6.0.tsv`
- `supporting/variant_changes_excluded_v1.6.0.tsv.gz`

## Files renamed

Every carried file is renamed `_v1.5.1` to `_v1.6.0`, its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.6.0.xlsx`. `supporting/site_changes_v1.5.1_to_v1.6.0.tsv` replaces v1.5.1's site changes table; it lists every one of the 11,632 methylarginine sites as kept.

Code: `scripts/66_stage_cancer_changes.py` (the K, S, T and Y change tables and the background), `scripts/67_stage_clinvar.py` (ClinVar), `scripts/69_stage_gnomad.py` (gnomAD), `scripts/68_release_cancer_changes.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at commit `36052256146a2c40b83c0e6ef2d58d11ae44f240`.
