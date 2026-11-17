# PTM_Pipeline v1.5.1: site-level interaction evidence for the context marks

Release v1.5.1 adds, for the site sets methyllysine, lysine, phospho, the site-level interaction tables that the methylarginine release has. Everything else is release v1.5.0, unchanged and not repeated here. The interaction tables of methylarginine are here too, under this release's names: they are the files of v1.5.0, byte for byte (which are those of v1.4.3), so one release holds the tables of every site set.

The tables are made by `scripts/64_export_site_interactions.py`, which runs the interface, crosslink and structure-contact steps of the methylarginine pipeline (`scripts/37_export_interaction_data.py`) with the modified residue as a parameter (`config/site_sets.json`). Run on the methylarginines it gives the methylarginine release tables (`tests/test_site_interactions.py`).

## Site sets

| Site set | Residues | Sites from | Sites | with a PDB interface | with RNA or DNA within 4 Å in a structure | with evidence from cells |
|---|---|---|---|---|---|---|
| methylarginine | R | `releases/v1.4.3/human_methylarginine_sites_v1.4.3.tsv` | 11,632 | 1,396 | 344 | 4,330 |
| methyllysine | K | `releases/v1.5.0/human_methyllysine_sites_v1.5.0.tsv` | 5,178 | 419 | 86 | 114 |
| lysine | K | `releases/v1.5.0/human_acetylation_sites_v1.5.0.tsv`; `releases/v1.5.0/human_lysine_acylation_sites_v1.5.0.tsv`; `releases/v1.5.0/human_ubiquitination_sites_v1.5.0.tsv`; `releases/v1.5.0/human_sumoylation_sites_v1.5.0.tsv` | 90,876 | 9,529 | 1,706 | 2,429 |
| phospho | S, T, Y | `releases/v1.5.0/human_phosphorylation_sites_v1.5.0.tsv` | 142,415 | 6,635 | 585 | 5,042 |

## Tables

One file per site set and table, `supporting/<site set>_<table>_v1.5.1.tsv`, in the columns of the methylarginine tables of the same name (`releases/v1.5.0/data_dictionary.tsv`). A site is `<accession>:<residue><position>`.

| Table | One row per |
|---|---|
| `site_interactions` | piece of evidence at a site: a PDB interface with a protein, RNA, DNA or another copy of the protein (PDBe-KB), an RNA or DNA crosslink, a CAPRI peptide, an RBR-ID residue, a UniProt DNA-binding region |
| `site_interaction_records` | site and record: each PDB entry or study behind `site_interactions`, with its method, paper and what the entry models at the site |
| `nucleic_acid_contacts` | site and RNA or DNA chain within 4 Å of its side chain in a PDB entry: the contacted nucleotides and the bound sequence around them |
| `nucleic_acid_proximal_sites` | site and complex with RNA or DNA within 8 Å of the side chain and none within 4 Å |
| `nucleic_acid_near_sites` | site and complex in which a residue at i ± 1–5, not the site, is within 4 Å of RNA or DNA (the side chain not within 8 Å) |
| `residue_evidence` | run of residues of a protein of the set with the same evidence |

Side-chain distances are measured from the atoms that carry the modification and the carbons next to them: lysine CD, CE, NZ; serine CB, OG; threonine CB, OG1, CG2; tyrosine CE1, CE2, CZ, OH.

### record_class

| Code | Meaning |
|---|---|
| `binds_modified` | binds with the site modified in the experiment: a PDB structure modelling a modified residue there (the components of config/site_sets.json modified_components) |
| `binds_unmodified` | binds with the unmodified residue at the site: a PDB structure modelling it there |
| `binds_other_residue` | binds with another residue at the site: a PDB structure of a mutated or variant construct (state_detail names it) |
| `binds_state_not_read` | a PDB interface at the site whose residue there could not be read (no entity of the entry aligns to the site) |
| `nucleic_acid` | RNA or DNA: contact in a PDB structure, a crosslink or RNA-bound peptide in cells, an RBR-ID residue, or a UniProt DNA-binding region (annotation) |
| `self` | another copy of the same protein in a PDB structure |

### site_state

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

The tables of methylarginine hold more kinds of evidence, classes and states than these (predicted interfaces and contacts, IntAct features, RBDmap and pCLAP): `releases/v1.5.0/data_dictionary.tsv` and its README describe them. In the table above, their evidence from cells counts those kinds too.

## Not in the tables of methyllysine, lysine, phospho

- RBDmap fragments and pCLAP peptides: their extents are not republished.
- Predicted interfaces and contacts (ProtVar, AlphaFold DB) and IntAct features: the pipeline computes them for methylarginines only.
- Enzyme-substrate records (kinases, writers and erasers).
- Statistics: no test is run on these tables.

## Sources and licences

| Source | Gives | Licence |
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

## Files

| File | Rows |
|---|---|
| `supporting/lysine_nucleic_acid_contacts_v1.5.1.tsv` | 28,471 |
| `supporting/lysine_nucleic_acid_near_sites_v1.5.1.tsv` | 7,652 |
| `supporting/lysine_nucleic_acid_proximal_sites_v1.5.1.tsv` | 10,820 |
| `supporting/lysine_residue_evidence_v1.5.1.tsv` | 179,341 |
| `supporting/lysine_site_interaction_records_v1.5.1.tsv` | 118,646 |
| `supporting/lysine_site_interactions_v1.5.1.tsv` | 19,011 |
| `supporting/methylarginine_nucleic_acid_contacts_v1.5.1.tsv` | 14,212 |
| `supporting/methylarginine_nucleic_acid_near_sites_v1.5.1.tsv` | 943 |
| `supporting/methylarginine_nucleic_acid_proximal_sites_v1.5.1.tsv` | 1,519 |
| `supporting/methylarginine_residue_evidence_v1.5.1.tsv` | 222,880 |
| `supporting/methylarginine_site_interaction_records_v1.5.1.tsv` | 52,060 |
| `supporting/methylarginine_site_interactions_v1.5.1.tsv` | 22,899 |
| `supporting/methyllysine_nucleic_acid_contacts_v1.5.1.tsv` | 1,466 |
| `supporting/methyllysine_nucleic_acid_near_sites_v1.5.1.tsv` | 678 |
| `supporting/methyllysine_nucleic_acid_proximal_sites_v1.5.1.tsv` | 1,326 |
| `supporting/methyllysine_residue_evidence_v1.5.1.tsv` | 58,802 |
| `supporting/methyllysine_site_interaction_records_v1.5.1.tsv` | 8,367 |
| `supporting/methyllysine_site_interactions_v1.5.1.tsv` | 1,358 |
| `supporting/phospho_nucleic_acid_contacts_v1.5.1.tsv` | 14,030 |
| `supporting/phospho_nucleic_acid_near_sites_v1.5.1.tsv` | 4,093 |
| `supporting/phospho_nucleic_acid_proximal_sites_v1.5.1.tsv` | 4,954 |
| `supporting/phospho_residue_evidence_v1.5.1.tsv` | 198,172 |
| `supporting/phospho_site_interaction_records_v1.5.1.tsv` | 81,420 |
| `supporting/phospho_site_interactions_v1.5.1.tsv` | 18,389 |
