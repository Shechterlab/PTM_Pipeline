# Human methylarginine sites, v2.0

The corrected union of human protein arginine methylation sites behind DeAngelo, Silverstein & Shechter, Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with the provenance of every site and of every source record. It replaces v1.0, the published union (PTM_results/, tag v1.0.0: 13,719 sites on 5,090 accessions), after the remapping corrections described in DATA_PROVENANCE_REVIEW.md.

- 13,353 sites on 4,720 reviewed human proteins, plus 10 sites included pending a policy decision (flag pending_splice_junction).
- 28,960 source records from 13 studies and databases, 27,968 of them placed.
- Since v1.0: 12,492 sites kept, 1,227 removed, 871 added.

## Files

| File | One row per | Use it for |
|---|---|---|
| human_methylarginine_sites_v2.0.tsv | site | the site list: sequence, methylation evidence, context, provenance |
| human_methylarginine_sites_v2.0.xlsx | (workbook) | all tables, with column notes, links and highlighting |
| human_methylarginine_source_records_v2.0.tsv | source record (one study's report of one site) | what each study reported and where it was placed |
| human_methylarginine_changes_v1.0_to_v2.0.tsv | site of v1.0 or v2.0 | what changed and why |
| human_methylarginine_studies_v2.0.tsv | study or database | citations, DOIs and how much each contributes |
| MD5SUMS | file | checksums (`md5sum -c MD5SUMS`) |

## Where the sites come from

Three source tables: Maron et al. 2021 Table S5 (a compilation of ten studies: Larsen 2016, PhosphoSitePlus, Musiani 2019, Guo 2014, Fong 2019, Maron 2021's own data, Li 2021, Fedoriw 2019, Lim 2020, Wei 2020), ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of the PTM master table (UniProt and IEDB records). Each record's peptide, window or position was placed on the reviewed canonical UniProt sequence by scripts/02_remap_arg_methyl_to_canonical.py and checked against its own evidence by scripts/qc_remap_evidence_check.py. Maron 2021 Tables S2 and S3 add methylation-degree calls, the antibodies that recovered each site and its response to PRMT inhibitors; they annotate sites but add none. A Table S2 or S3 site is attached only when its peptide sits on the reference entry at the reported start, or occurs once in the entry (2,880 of 2,884 Table S2 rows, 24,680 of 24,696 Table S3 site rows).

## Site table columns

**Site**

- `site_id`: UniProt accession and arginine position, e.g. P35637:R216.
- `gene`: primary gene name of the UniProt entry.
- `protein_name`: recommended protein name of the UniProt entry.
- `uniprot_accession`: reviewed canonical UniProtKB accession (linked to its UniProt page in the workbook).
- `uniprot_entry_name`: UniProtKB entry name.
- `position`: position of the arginine in the reference sequence (1-based).
- `residue`: residue at that position in the reference sequence (always R).

**Sequence**

- `window_15`: residues -7 to +7 around the site; _ beyond a protein end.
- `window_31`: residues -15 to +15 around the site; _ beyond a protein end.
- `protein_length`: length of the reference sequence.
- `reference_sequence_md5`: MD5 of the reference sequence (UniProt sequence MD5 format: uppercase hex); compare with a later UniProt entry to see whether the numbering still holds.

**Methylation evidence**

- `methyl_states_reported`: methylation degree given by the source records: me1 (monomethyl), me2 (dimethyl, asymmetric or symmetric not stated); blank where the sources give none.
- `maron2021_ip_states`: degree called in Maron 2021 Table S2 (PTMScan methyl-arginine immunoprecipitation): me1, me2.
- `maron2021_ip_antibodies`: antibodies whose immunoprecipitation recovered the site in Maron 2021 Table S3: ADMA (asymmetric dimethyl), SDMA (symmetric dimethyl), MMA (monomethyl).
- `maron2021_GSK591_PRMT5i`: response to the type II (PRMT5) inhibitor GSK591 in Maron 2021 Table S3: decreased or increased (adjusted p < 0.05 in at least one immunoprecipitation), mixed, or no significant change.
- `maron2021_MS023_typeIi`: response to the type I PRMT inhibitor MS023 in Maron 2021 Table S3, as for GSK591.
- `maron2021_inhibitor_detail`: per inhibitor and immunoprecipitation: log2 intensity difference vs DMSO and adjusted p value of the strongest peptide (lowest adjusted p). Measured on peptides: arginines on the same peptide share its values.

**Context**

- `motifs`: sequence motifs around the site, as defined by scripts/08: RG, RGG, GRG (G at -1 and +1), GAR (G at -2, A at -1), RXR (R at -2 or +2), CARM1_proline_rich (P within -5..+5), CARM1_hydrophobic_proline (P within -5..+5 and a hydrophobic residue within -3..+3), DR, ER (D or E at -1), RD, RE (D or E at +1).
- `disorder`: MobiDB predicted disorder at the site (scripts/06): disordered; disorder boundary (outside, within 20 residues of a disordered region); ordered; no disorder annotation.
- `disorder_edge_distance`: residues to the nearest edge of a predicted disordered region.
- `domain`: InterPro domain or repeat context (scripts/05): in domain; boundary (within 20 residues of one); inter-domain linker; distal; no domain annotation; not staged (protein outside the InterPro intervals staged for the published union).
- `nearest_domain`: InterPro entry nearest the site (the innermost when inside several), with its accession and span.
- `domain_edge_distance`: residues to the nearest domain edge (0 inside a domain).

**Provenance**

- `n_studies`: number of studies or databases reporting the site.
- `studies`: the studies or databases (one row each in the Studies table).
- `study_dois`: their DOIs, in the same order.
- `n_records`: number of source records placed on the site.
- `n_records_at_reported_position`: records whose own accession and position name this site unchanged.
- `source_tables`: source tables of the records: Maron 2021 Table S5, ProMetheusDB mmc2, PTM master table.
- `record_ids`: the supporting rows of the source records table (record_id).

**Status**

- `release_status`: placed, or included pending a policy decision (see flags).
- `flags`: cautions, separated by ";": pending_splice_junction; paralog_window; near_match_only; shared_sequence; single_study.
- `same_window_in`: other reviewed entries holding the same 31 residues (see the shared_sequence flag).

**Change**

- `in_previous_version`: whether the previous version (the published union) had this site.
- `change_since_previous_version`: kept or added; removed sites are listed in the changes table.
- `change_reasons`: for an added site, why (reason codes, see the Reason codes table).

## Flags

- `pending_splice_junction` (10 sites): a window crossing a splice junction: one full flank (at least 6 residues) matches the entry exactly at the reported arginine, the other does not, and no other entry holds the window. The pipeline leaves it unplaced; the release includes it until the policy is decided.
- `paralog_window` (22 sites): a supporting window matches a paralog exactly and this entry within 2 residues (issue 9, open; accepted in config/remap_evidence_allowlist.tsv).
- `near_match_only` (12 sites): every supporting record matches the entry within 2 residues, not exactly (a variant or sequence conflict).
- `shared_sequence` (593 sites): residues -15 to +15 occur identically in another reviewed entry (listed in same_window_in): histone variants, identical paralogs. The evidence cannot tell the proteins apart; the gene is the one the source named.
- `single_study` (6,450 sites): one study or database reports the site.

Methylation degree is known for 2,473 sites (sources or Maron Table S2); 957 sites have Maron 2021 antibody and inhibitor data. Domain context is "not staged" for 126 sites on 94 proteins that were not in the published union, because the frozen InterPro intervals were staged for that union only.

## Source records, changes and studies

The columns of the other tables are listed on the workbook's Columns sheet. Every record carries the code of the reason its placement changed; every removed or added site carries the codes of its records:

| reason_code | issue | records | meaning |
|---|---|---|---|
| peptide_recentred | 1 | 1,773 | the published run moved the peptide to the arginine nearest its middle; the release keeps the reported arginine inside the peptide |
| unreviewed_accession | 2 | 879 | the published run filed reviewed-entry numbering under the source's unreviewed or obsolete accession; the release files it under the reviewed entry |
| lysine_window | 3 | 132 | the source window is centred on a lysine; the published run turned it into a nearby arginine site |
| truncated_window | 4 | 303 | the window is cut short at a protein end; the published run read it from the wrong end |
| multi_arginine_peptide | 5 | 54 | a peptide with several arginines that does not cover the reported position; the published run guessed one |
| wrong_protein | 6 | 28 | the published run filed the position under another reviewed protein than the one whose sequence matched |
| outdated_accession_or_symbol | 7 | 103 | the source accession is merged or obsolete, or its gene symbol outdated; the release reaches the current entry |
| evidence_mismatch | 8 | 88 | the published run kept the reported arginine although the source peptide or window does not match the entry there |
| window_padding | 13 | 5 | the window is padded where it runs off a protein end; the padding now gives the cut side |
| ambiguous_evidence | 14 | 0 | the evidence supports several sites equally; the release leaves the record unplaced |
| other |  | 0 | placement changed for a reason not covered by the codes above |

## How it was built

- Code: PTM_Pipeline commit 76a3d8e4ca473c4e1dc90d1fdcfa23aed0c29ed1; driver hpc/build_release.sh (parse, place, union, the evidence check and the residue audit, which stop the build on any failure), then scripts/build_release.py.
- Inputs: inputs/, checked against inputs/MANIFEST.tsv. Reference sequences: UniProtKB reviewed human (proteome UP000005640), staged 2026-04-07; release not recorded, earlier than 2026_03. The published union was mapped to the same sequences. Accession history: UniProt 2026_03.
- Context: MobiDB disorder (export of 2026-03-30), InterPro domains and repeats (staged 2026-04-07), motifs, with the classifiers of scripts/06, scripts/05 and scripts/08.

## Licences

PhosphoSitePlus® data (www.phosphosite.org) are licensed under CC BY-NC-SA 3.0: non-commercial use, with attribution, adaptations shared under the same licence. The release as a whole is CC BY-NC-SA 4.0 (LICENSE.md). Cite the studies in study_dois for the sites you use, and the paper above for the union.

## Run record

```
built: 2026-10-20T01:24:40Z
commit: 76a3d8e4ca473c4e1dc90d1fdcfa23aed0c29ed1
uncommitted changes:
python: 3.11.17 pandas 2.2.3 numpy 2.2.6 openpyxl 3.1.5
inputs: inputs (all match MANIFEST.tsv)
```
