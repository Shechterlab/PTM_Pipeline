# Human methylarginine sites: release 2026-10-corrected

The corrected union of human protein arginine methylation sites behind DeAngelo, Silverstein & Shechter, Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with the provenance of every site and of every source record. It replaces the published union (PTM_results/, 13,719 sites on 5,090 accessions) after the remapping corrections described in DATA_PROVENANCE_REVIEW.md; changes_from_published.tsv lists every difference.

## Contents

- 13,353 sites on 4,720 reviewed human proteins, placed on the reference sequences below; plus 10 sites included pending a policy decision (flag pending_splice_junction).
- 28,960 source records, 27,968 of them placed.
- Against the published union: 12,492 sites kept, 1,227 removed, 871 added.

## Files

| File | One row per | Use it for |
|---|---|---|
| sites.tsv | site | the site list, with where each site comes from |
| evidence.tsv | source record (a study's report of a site) | what each study reported and where it was placed |
| changes_from_published.tsv | site of the published union or of this release | what changed and why |
| methylarginine_sites.xlsx | the three tables and this README | the same, as one workbook |
| release_manifest.tsv | file | size and MD5 of each file |

## Where the sites come from

Three source tables: Maron et al. 2021 Table S5 (a compilation of ten studies: Larsen 2016, PhosphoSitePlus, Musiani 2019, Guo 2014, Fong 2019, Maron 2021's own data, Li 2021, Fedoriw 2019, Lim 2020, Wei 2020), ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of the PTM master table (UniProt and IEDB records). Each record's peptide, window or position was placed on the reviewed canonical UniProt sequence by scripts/02_remap_arg_methyl_to_canonical.py and checked by scripts/qc_remap_evidence_check.py.

## sites.tsv columns

- site, accession, entry_name, gene, protein_name, position, residue: the site on the reviewed canonical entry.
- window_31: residues -15 to +15 of the reference sequence (_ beyond the ends).
- protein_length, reference_sequence_md5: the reference sequence the position refers to (MD5 in UniProt's sequence.md5 format, to compare with the current UniProt entry).
- release_status: placed, or pending a policy decision (see flags).
- flags: single_study (one study reports the site); paralog_window (a supporting window matches a paralog exactly and this entry within 2 residues; issue 9, open); near_match_only (every supporting record matches the entry within 2 residues, not exactly); pending_splice_junction (a window crossing a splice junction: one full flank matches exactly at the reported arginine, the other does not; the pipeline leaves it unplaced, the release includes it until the policy is decided).
- methyl_state, methyl_states_reported: methylation state calls, where a source gives them.
- n_records, n_studies, studies, study_dois, source_tables: the support of the site.
- n_records_at_reported_position: records whose own accession and position name this site unchanged.
- in_published_union: whether the published union had this site.
- record_ids: the supporting rows of evidence.tsv.

## evidence.tsv columns

- record_id; source_table; study, study_doi, reference_text, pmid: who reported it.
- reported_accession, reported_gene, reported_position, reported_methyl_state: what the source reported.
- evidence_kind, evidence_sequence: the window, peptide or marked peptide the source gave, verbatim.
- placed_accession, placed_position, placed_site, placement, placement_code: where the release puts it and how.
- evidence_check: the independent check (agrees; contradicts evidence; one of several evidence sites; ...).
- accepted_exception: why a contradiction is accepted for now (open issue 9).
- published_site, change_from_published, reason_code, issue: where the published run put it and why it differs.

## Reason codes

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

Flags in sites.tsv: single_study 6,450, paralog_window 22, near_match_only 12, pending_splice_junction 10.

## How it was built

- Code: PTM_Pipeline commit 9bb6716065038467e2ce799f15235ffa006a9a9c (scripts/02, scripts/qc_remap_evidence_check.py, scripts/build_release.py; driver hpc/build_release.sh).
- Inputs: inputs/ (checksums in inputs/MANIFEST.tsv). Reference sequences: UniProtKB reviewed human (proteome UP000005640), staged 2026-04-07; release not recorded, earlier than 2026_03. The published union was mapped to the same sequences. Accession history: UniProt 2026_03.
- Checks passed: every record agrees with its own evidence except the accepted issue-9 rows in config/remap_evidence_allowlist.tsv; every placed site is an arginine of the reference sequence.

## Licences

PhosphoSitePlus® data (www.phosphosite.org) are licensed under CC BY-NC-SA 3.0: non-commercial use, with attribution, adaptations shared under the same licence. The release as a whole is CC BY-NC-SA 4.0 (LICENSE.md). Cite the studies in study_dois for the sites you use.

## Run record

```
built: 2026-10-20T00:22:35Z
commit: 9bb6716065038467e2ce799f15235ffa006a9a9c
uncommitted changes:
python: 3.11.17 pandas 2.2.3 numpy 2.2.6
inputs: inputs (all match MANIFEST.tsv)
```
