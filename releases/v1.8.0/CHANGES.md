# Changes from v1.7.0 to v1.8.0

Everything of v1.7.0 is in v1.8.0 under this release's names. The owner's decisions of 23 November 2026: the per-site DepMap cell lines and AlphaMissense scores go in separate per-site tables (the site tables keep their bytes); ClinVar gets disease families; the site-set table gets the scope of its test; the title names every mark. AlphaMissense is a new source, so this is a minor release.

## New: per-site annotations, ClinVar disease families, and the scope of the site-set test

- **`supporting/site_annotations_<set>_v1.8.0.tsv`** (methylarginine, methyllysine, lysine, phospho): one row per site of the set (accession, position, residue). Each row gives the DepMap 24Q4 cell lines with a missense or stop change at the site (count, names and model IDs), and AlphaMissense's mean and highest predicted pathogenicity over the residue's 19 substitutions, with how many of them it classes pathogenic. The site tables themselves keep their bytes (the methylarginine table is unchanged since v1.4.3, the context-mark tables since v1.5.0); join the annotations to them on accession, position and residue.
- **AlphaMissense** counts for an entry when its sequence equals the 2026_03 canonical sequence over its whole length: 19,382 entries do, 851 differ and 183 are missing (`inputs/layers/alphamissense_status.tsv`). AlphaMissense (Cheng et al., Science 381:eadg7492, 2023; predictions doi:10.5281/zenodo.10813168, CC BY 4.0 by the Zenodo record, the file header still naming CC BY-NC-SA 4.0) predicts the pathogenicity of each amino-acid substitution; its scores are less reliable in intrinsically disordered regions, where many modification sites lie. 0.564 is AlphaMissense's threshold for calling a single substitution likely pathogenic (Cheng et al. 2023); a residue whose mean over all 19 lies above it is one where most substitutions are predicted damaging. Cell-line names and model IDs are separated by | (one name, RS4;11, holds a semicolon).
- **`disease_families`** is appended to `supporting/clinvar_changes_<residue>_v1.8.0.tsv.gz`. Each variation's conditions are classed by `classify_clinvar()` in `scripts/disease_families.py`: ClinVar-specific rules first, then the DisGeNET rules of `classify()`, which is unchanged (so the disease analyses of earlier releases stay reproducible). On ClinVar names `classify()` alone misfiles common conditions (Lynch and Li-Fraumeni syndromes as developmental, retinal dystrophies as neurodegeneration, the developmental and epileptic encephalopathies as craniofacial, malignant hyperthermia as cancer); the ClinVar rules correct these and add a family, Epilepsy and seizures. Of the 43,394 pathogenic or likely pathogenic rows with a named condition, 36,242 (83.5%) get a family.
- **`sites_total`, `sites_weighted`, `k_total`** are appended to `supporting/site_set_propensity_v1.8.0.tsv`: the test covers sites_weighted of the sites_total sites and observed of the k_total changes (the rest lie on proteins without weights). Every earlier column keeps its values.

| Site set | Sites | with a DepMap cell line changed | with AlphaMissense | mean of the 19 substitutions above 0.564 |
|---|---:|---:|---:|---:|
| Methylarginine (R) | 11,632 | 849 | 11,273 | 5,904 |
| Methyllysine (K) | 5,178 | 151 | 4,894 | 2,313 |
| Lysine acylation, ubiquitination and SUMOylation (K) | 90,876 | 2,466 | 87,422 | 43,840 |
| Phosphorylation (S, T, Y) | 142,415 | 5,303 | 137,504 | 35,412 |

The title of the release now names every mark it holds.

## New files

- `supporting/site_annotations_lysine_v1.8.0.tsv`
- `supporting/site_annotations_methylarginine_v1.8.0.tsv`
- `supporting/site_annotations_methyllysine_v1.8.0.tsv`
- `supporting/site_annotations_phospho_v1.8.0.tsv`

## Columns appended

- `supporting/clinvar_changes_K_v1.8.0.tsv.gz`: `disease_families`
- `supporting/clinvar_changes_R_v1.8.0.tsv.gz`: `disease_families`
- `supporting/clinvar_changes_S_v1.8.0.tsv.gz`: `disease_families`
- `supporting/clinvar_changes_T_v1.8.0.tsv.gz`: `disease_families`
- `supporting/clinvar_changes_Y_v1.8.0.tsv.gz`: `disease_families`
- `supporting/site_set_propensity_v1.8.0.tsv`: `sites_total`, `sites_weighted`, `k_total`

## Files renamed

Every carried file is renamed `_v1.7.0` to `_v1.8.0`, its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.8.0.xlsx`, and the tables above with columns appended. `supporting/site_changes_v1.7.0_to_v1.8.0.tsv` replaces v1.7.0's site changes table; it lists every one of the 11,632 methylarginine sites as kept.

Code: `scripts/72_stage_alphamissense.py` (AlphaMissense per residue), `scripts/disease_families.py` (`classify_clinvar`), `scripts/propensity.py` (`set_scope`), `scripts/73_release_site_annotations.py` (this release), `scripts/65_carry_release.py`, at commit `461054cb44aa22e6902d351564f039d6015eaf4c`.
