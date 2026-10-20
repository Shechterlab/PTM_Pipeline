# Changes from v1.0.0 to v1.1.0

v1.1.0 corrects how source records were placed on the reference protein sequences. The sources, the reference sequences and the analysis are those of v1.0.0, the published union (repository tag v1.0.0).

| | v1.0.0 | v1.1.0 |
|---|---|---|
| Sites | 13,719 | 13,358 (and 10 pending, flagged) |
| Proteins | 5,090 | 4,722 |
| Source records placed | 28,134 of 28,960 | 27,983 of 28,960 |

12,495 sites are kept, 1,224 removed and 873 added. A site moves when its records now point at another arginine: the old one is removed and the new one added.

## What was corrected

Each record whose placement changed carries one reason code. Issue numbers refer to DATA_PROVENANCE_REVIEW.md.

| Reason code | Issue | What changed | Records | Sites removed | Sites added |
|---|---|---|---|---|---|
| peptide_recentred | 1 | the published run moved the peptide to the arginine nearest its middle; the release keeps the reported arginine inside the peptide | 1,773 | 117 | 482 |
| unreviewed_accession | 2 | the published run filed reviewed-entry numbering under the source's unreviewed or obsolete accession; the release files it under the reviewed entry | 879 | 857 | 268 |
| lysine_window | 3 | the source window is centred on a lysine; the published run turned it into a nearby arginine site | 132 | 71 | 0 |
| truncated_window | 4 | the window is cut short at a protein end; the published run read it from the wrong end | 303 | 138 | 92 |
| multi_arginine_peptide | 5 | a peptide with several arginines that does not cover the reported position; the published run guessed one | 39 | 5 | 0 |
| peptide_start | 5 | a peptide with several arginines, numbered differently from the entry; the published run took the arginine nearest its middle; this version counts the reported position from the peptide's start in Maron 2021 Tables S2/S3 | 8 | 1 | 2 |
| wrong_protein | 6 | the published run filed the position under another reviewed protein than the one whose sequence matched | 28 | 23 | 13 |
| outdated_accession_or_symbol | 7 | the source accession is merged or obsolete, or its gene symbol outdated; the release reaches the current entry | 103 | 0 | 30 |
| evidence_mismatch | 8 | the published run kept the reported arginine although the source peptide or window does not match the entry there | 88 | 21 | 0 |
| window_padding | 13 | the window is padded where it runs off a protein end; the padding now gives the cut side | 5 | 0 | 1 |

A site can carry several codes when its records changed for different reasons.

## Proteins featured in the paper

| Protein | Sites in v1.0.0 | Sites in v1.1.0 | Removed | Added |
|---|---|---|---|---|
| FUS (P35637) | 30 | 34 | - | R422, R518, R522, R524 |
| EWSR1 (Q01844) | 35 | 36 | - | R518 |
| FMR1 (Q06787) | 13 | 13 | R562, R564 | R483, R558 |
| Histone H3.1 (P68431) | 10 | 10 | - | - |
| SmB/SNRPB (P14678) | 13 | 13 | - | - |
| SmD1/SNRPD1 (P62314) | 6 | 14 | - | R92, R104, R108, R110, R112, R114, R118, R119 |
| SmD3/SNRPD3 (P62318) | 6 | 7 | - | R125 |

## Decided in this version

- **Peptides with several arginines.** A peptide does not say which of its arginines is methylated, and the position the source reports may follow an older numbering. Where Maron et al. 2021 Tables S2 and S3 give the peptide's start in the same numbering, the reported position is counted from that start (reason peptide_start): 15 records, on ARID1B R608, ARID1B R640, ARID1B R648, CENPVL3 R17, TENT4A R145, WWC3 R493, WWC3 R500, WWC3 R507. Without such a start the record stays unplaced (reason multi_arginine_peptide). Where a start is known, counting from it agrees with all 2,949 placements the pipeline made from plain peptides by other means.
- **Domain context for every protein.** The published analyses used InterPro intervals staged for only part of the published union: 573 of its reviewed proteins had none, so their 847 sites were classed "no domain annotation", although InterPro has a domain or repeat on 457 of these proteins (681 sites). This version stages InterPro 108.0, the release current when the reference sequences were staged, for all 20,416 reviewed entries; for the proteins the published file covers it gives the same intervals. The domain context of 658 kept sites changes as a result.

## Still open

- **Windows crossing a splice junction** (10 sites, flag pending_splice_junction; 10 of them in v1.0.0): one full flank of the window matches the entry at the reported arginine and the other does not, and no other entry holds the window. The pipeline leaves these records unplaced; the release keeps the sites, flagged, until the policy is decided: CTTNBP2 R451, DIS3 R144, DIS3 R147, DIS3 R150, GANAB R185, GCC2 R303, RUSC2 R882, SNX12 R153, SQSTM1 R312, TCF7L2 R126.
- **Paralog windows** (22 sites, flag paralog_window; issue 9): a window matches a paralog exactly and the named entry within 2 residues. The records stay on the entry the source named, accepted in config/remap_evidence_allowlist.tsv.

## Tracing a change

- supporting/site_changes_v1.0.0_to_v1.1.0.tsv: one row per site of either version, with its change, reason codes, the records involved and, for a removed site, where those records are now (for an added site, where they were).
- supporting/source_records_v1.1.0.tsv: one row per source record, with what the study reported, the evidence sequence, where each version placed it, the independent evidence check and the reason code.

Columns of the source records table:

- `record_id`: row number of the record, cited by record_ids in the other tables.
- `source_table`: the source table the record comes from.
- `study`: the study or database that reported it.
- `study_doi`: its DOI.
- `reference_text`: the reference as the source table gives it.
- `pmid`: PubMed ID, where the source gives one.
- `reported_accession`: accession as reported.
- `reported_gene`: gene as reported.
- `reported_position`: position as reported.
- `reported_methyl_state`: methylation state as reported.
- `evidence_kind`: kind of sequence evidence: centred window (length), peptide, marked peptide, or position only.
- `evidence_sequence`: the window, peptide or marked peptide, verbatim.
- `placed_site`: where this version places the record (blank: not placed).
- `placed_accession`: accession of the placed site.
- `placed_position`: position of the placed site.
- `placement`: how the record was placed, in words.
- `placement_code`: the same, as scripts/02 remap_status.
- `evidence_check`: independent check of the placement against the record's own evidence (scripts/qc_remap_evidence_check.py).
- `accepted_exception`: why a contradiction is accepted for now (open issue 9).
- `previous_site`: where the previous version placed the record.
- `change_since_previous_version`: same, moved, newly placed, no longer placed, or not placed in either.
- `reason_code`: why the placement changed (CHANGES.md).
- `issue`: the issue number of DATA_PROVENANCE_REVIEW.md behind the reason.

Columns of the site changes table:

- `site_id`: the site.
- `gene`: primary gene name.
- `uniprot_accession`: reviewed canonical accession.
- `position`: arginine position.
- `in_previous_version`: site of the previous version.
- `in_this_version`: site of this version.
- `change`: kept, removed or added.
- `reason_codes`: why it was removed or added (CHANGES.md).
- `issues`: issue numbers of DATA_PROVENANCE_REVIEW.md.
- `records_now_at_or_previously_at`: for a removed site, where its records are placed now; for an added site, where they were placed before.
- `record_ids`: the records involved (source records table).
