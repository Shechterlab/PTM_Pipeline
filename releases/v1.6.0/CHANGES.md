# Changes from v1.5.1 to v1.6.0

Everything of v1.5.1 is in v1.6.0 under this release's names, its data unchanged. v1.6.0 adds the somatic changes at every candidate residue of the site sets, in TCGA tumors and DepMap cell lines.

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

## New files

- `supporting/cancer_change_protein_background_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_change_models_v1.6.0.tsv`
- `supporting/depmap_24q4_changes_K_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_R_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_S_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_T_v1.6.0.tsv.gz`
- `supporting/depmap_24q4_changes_Y_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_K_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_R_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_S_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_T_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_changes_Y_v1.6.0.tsv.gz`
- `supporting/tcga_mc3_projects_v1.6.0.tsv`

## Files renamed

Every carried file is renamed `_v1.5.1` to `_v1.6.0`, its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.6.0.xlsx`. `supporting/site_changes_v1.5.1_to_v1.6.0.tsv` replaces v1.5.1's site changes table; it lists every one of the 11,632 methylarginine sites as kept.

Code: `scripts/66_stage_cancer_changes.py` (the K, S, T and Y change tables and the background), `scripts/68_release_cancer_changes.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at commit `72818af85bf146475d676a7af35c362cf62fea78`.
