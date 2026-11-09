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

Versioned tables of other marks, from openly licensed sources, for sites of those marks to read (`README.md` has the sources, licenses, counts and how to cite).

- `human_methyllysine_sites_v1.5.0.tsv`: 756 rows on 412 residues of 190 proteins; by state me 251, me1 256, me2 127, me3 122; by evidence class by-similarity 190, experimental 360, unspecified 206; residues by source UniProt only 178, both 219, iPTMnet only 15.
- `human_phosphorylation_sites_v1.5.0.tsv`: 72,134 rows on 72,134 residues of 10,525 proteins; by state phospho 72,134; by evidence class by-similarity 7,802, experimental 61,719, text-mined 2,526, unspecified 87; residues by source UniProt only 8,017, both 32,911, iPTMnet only 31,206.
- `human_methylarginine_crosscheck_v1.5.0.tsv`: 2,471 arginine methylation rows of UniProt and iPTMnet with whether the methylarginine data has the arginine; residues by source UniProt only 740, both 742, iPTMnet only 65.

Sources: UniProt (release 2026_04, CC BY 4.0) and iPTMnet (release 6.2, CC BY-NC-SA 4.0). The release stays CC BY-NC-SA 4.0. iPTMnet's download readme also mentions CC BY 4.0; the stricter terms are followed.

Code: `scripts/54_stage_uniprot_ptm_sites.py`, `scripts/55_stage_iptmnet_ptm_sites.py`, `scripts/56_build_context_marks.py`; staged tables in `inputs/marks/`.
