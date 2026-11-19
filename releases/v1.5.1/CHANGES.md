# Changes from v1.5.0 to v1.5.1

v1.5.1 is the release; v1.5.0 was prepared but not released on its own. Everything of v1.5.0 is in v1.5.1 under this release's names, its data unchanged, and v1.5.1 adds site-level interaction tables for the context marks.

## New: site-level interaction evidence for the context marks

For methyllysine, lysine, phospho, `supporting/<site set>_<table>_v1.5.1.tsv` gives what the methylarginine release gives for its sites: PDB interfaces with each entry's method, paper and the state of the residue it models, RNA and DNA contacts within 4 Å in structures with the bound sequences, proximal and nearby nucleic-acid contacts, and RNA-binding evidence from cells (`scripts/64_export_site_interactions.py`, the steps of `scripts/37` with the modified residue as a parameter).

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

## The methylarginine interaction tables, carried over

The six methylarginine interaction tables of v1.5.0 (and of v1.4.3) are carried over byte for byte under the names the other site sets use, so that every set's tables sit side by side:

| In v1.5.0 | In v1.5.1 |
|---|---|
| `supporting/site_interactions_v1.5.0.tsv` | `supporting/methylarginine_site_interactions_v1.5.1.tsv` |
| `supporting/residue_evidence_v1.5.0.tsv` | `supporting/methylarginine_residue_evidence_v1.5.1.tsv` |
| `supporting/nucleic_acid_contacts_v1.5.0.tsv` | `supporting/methylarginine_nucleic_acid_contacts_v1.5.1.tsv` |
| `supporting/nucleic_acid_near_sites_v1.5.0.tsv` | `supporting/methylarginine_nucleic_acid_near_sites_v1.5.1.tsv` |
| `supporting/nucleic_acid_proximal_sites_v1.5.0.tsv` | `supporting/methylarginine_nucleic_acid_proximal_sites_v1.5.1.tsv` |
| `supporting/site_interaction_records_v1.5.0.tsv` | `supporting/methylarginine_site_interaction_records_v1.5.1.tsv` |

## Files renamed

Every carried file is renamed `_v1.5.0` to `_v1.5.1`, its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.5.1.xlsx`. `supporting/site_changes_v1.5.0_to_v1.5.1.tsv` replaces v1.5.0's site changes table; it lists every one of the 11,632 methylarginine sites as kept.

## Before the release: CPLM 4.0 taken out of v1.5.0

CPLM 4.0 states no license ("All Rights Reserved") and compiles sites from other databases as well as from papers. The owner decided on 17 November 2026 to leave it out, and v1.5.0, not yet released, was rebuilt without it (pull request #27). Ochoa 2020, Berryhill 2023, Hendriks 2017 and Akimov 2018 stay. The context-mark tables of this release:

| Table | Rows with CPLM | Rows in this release |
|---|---:|---:|
| `human_methyllysine_sites_v1.5.1.tsv` | 10,889 | 5,746 |
| `human_lysine_acylation_sites_v1.5.1.tsv` | 88,816 | 8,008 |
| `human_acetylation_sites_v1.5.1.tsv` | 53,177 | 4,597 |
| `human_ubiquitination_sites_v1.5.1.tsv` | 125,771 | 63,480 |
| `human_sumoylation_sites_v1.5.1.tsv` | 48,593 | 40,249 |
| `human_phosphorylation_sites_v1.5.1.tsv` | 142,415 | 142,415 |

The site-level interaction tables of methyllysine and lysine were rebuilt from these tables; each holds exactly the rows of the sites still in the set.

## Fixed: the analyses checksums of v1.5.0

`analyses/MD5SUMS` and `analyses_with_lower_confidence/MD5SUMS` listed their README.md with the checksum of the v1.4.3 text, while the README names the release's own version; they are checksummed again here (and in v1.5.0).

Code: `scripts/64_export_site_interactions.py` (interaction tables), `scripts/56_build_context_marks.py` (context marks), `scripts/65_carry_release.py` (this release), at commit `854071a915e8a08f9090360c55ecb29efd96ae72`.
