# Changes from v1.4.3 to v1.5.0

v1.5.0 adds tables of context marks next to the methylarginine data. **The methylarginine data is unchanged from v1.4.3**: its 100 files are the v1.4.3 files under this release's names, with the data not rebuilt (`release.json`, `carried_over`, gives the v1.4.3 name and MD5 of each), because the build (`hpc/build_release.sh`) needs R with clusterProfiler and DOSE and the network, and no methylarginine input changed. Sites, proteins, records, statistics and layers are those of v1.4.3.

### Files renamed

Every file name carries its release's version, so each of the 100 carried files is renamed `_v1.4.3` to `_v1.5.0` (for example `human_methylarginine_sites_v1.4.3.tsv` is now `human_methylarginine_sites_v1.5.0.tsv`), its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.5.0.xlsx`. `supporting/site_changes_v1.4.3_to_v1.5.0.tsv` replaces v1.4.3's (v1.4.2 to v1.4.3); it lists every one of the 11,632 sites as kept, none added or removed, with its records the same.

| | v1.4.3 | v1.5.0 |
|---|---:|---:|
| Sites | 11,632 | 11,632 |
| Proteins | 4,439 | 4,439 |
| Sites kept pending | 0 | 0 |
| Source records | 67,516 | 67,516 |
| Studies and databases | 19 | 19 |

## New: context marks

Versioned tables of other marks, for sites of those marks to read (`README.md` has the sources, licenses, counts and how to cite). The sources are UniProt and iPTMnet, which are openly licensed, and, for methyllysine only, CPLM 4.0, which is not (below).

- `human_methyllysine_sites_v1.5.0.tsv`: 5,899 rows on 5,481 residues of 2,967 proteins; by state me 5,394, me1 256, me2 127, me3 122; by evidence class by-similarity 190, experimental 2,849, unspecified 2,860; residues by source CPLM only 5,069, UniProt + CPLM 75, UniProt + iPTMnet 8, UniProt + iPTMnet + CPLM 211, UniProt only 103, iPTMnet + CPLM 1, iPTMnet only 14.
- `human_phosphorylation_sites_v1.5.0.tsv`: 72,134 rows on 72,134 residues of 10,525 proteins; by state phospho 72,134; by evidence class by-similarity 7,802, experimental 61,719, text-mined 2,526, unspecified 87; residues by source UniProt + iPTMnet 32,911, UniProt only 8,017, iPTMnet only 31,206.
- `human_methylarginine_crosscheck_v1.5.0.tsv`: 2,471 arginine methylation rows of UniProt and iPTMnet with whether the methylarginine data has the arginine; residues by source UniProt + iPTMnet 742, UniProt only 740, iPTMnet only 65.

Sources: UniProt (release 2026_04, CC BY 4.0), iPTMnet (release 6.2, CC BY-NC-SA 4.0) and CPLM 4.0 (below). The release stays CC BY-NC-SA 4.0. iPTMnet's download readme also mentions CC BY 4.0; the stricter terms are followed.

### CPLM 4.0: license note

The methyllysine table includes all 5,356 validated human lysine methylation sites of CPLM 4.0 (2,939 proteins; Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849; `Methylation.txt` of 2021-06-13). **CPLM states no license**: cplm.biocuckoo.cn carries only a copyright notice ("All Rights Reserved"). The sites are included by the owner's decision of 2026-11-09 (confirmed 2026-11-10), cited, **pending the authors' permission**; the terms the authors give, when they answer, are followed. The CC BY-NC-SA 4.0 of the release does not grant rights in CPLM's records beyond what the authors allow. Each CPLM row has `CPLM` in `sources` and CPLM's own evidence label in `cplm_evidence` (Dat. 2,858, Exp. 618, Exp./Dat. 1,880), so a reader can filter on it or leave CPLM out. CPLM rows are state `me` (CPLM gives no degree); `evidence_class` is experimental for Exp. and Exp./Dat., and unspecified for Dat. (a record CPLM took from another database). `pmids` holds the references that are PubMed IDs (2,919 rows; the rest of CPLM's reference field is identifiers of other databases). Dropped on validation, nothing remapped: 494 accession not in canonical fasta, 60 canonical not K, 5 window mismatch. No lysine acylation is part of this release.

Code: `scripts/39_stage_cplm4_lys_methylation.py`, `scripts/54_stage_uniprot_ptm_sites.py`, `scripts/55_stage_iptmnet_ptm_sites.py`, `scripts/56_build_context_marks.py`; staged tables in `inputs/marks/`.
