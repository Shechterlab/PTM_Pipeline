# Changes from v1.0.0 to v1.1.0

v1.1.0 uses the sources and methods of DeAngelo et al. (2026), mapped to UniProt release 2026_03, with corrections to how some source records were placed on UniProt sequences. v1.0.0 is the published union (git tag v1.0.0, `releases/v1.0.0/`).

| | v1.0.0 | v1.1.0 |
|---|---:|---:|
| Sites | 13,719 | 13,358, and 10 pending |
| Proteins | 5,090 | 4,722 |
| Source records placed | 28,134 of 28,960 | 27,983 of 28,960 |

Of the 13,719 published sites, 12,495 are unchanged and 1,224 were removed; 873 sites were added. A site moves when its records now point at another arginine: the old one is removed and the new one added. One kept site has a new number only: RAB12 R10 is R106 in UniProt 2026_03, which adds 96 residues at its N-terminus.

## What was corrected

Each record whose placement changed carries one reason code.

| Reason code | What changed | Records | Sites removed | Sites added |
|---|---|---:|---:|---:|
| `peptide_recentred` | the previous version moved the peptide to the arginine nearest its middle. This version keeps the reported arginine inside the peptide | 1,773 | 117 | 482 |
| `unreviewed_accession` | the previous version filed reviewed-entry numbering under the source's unreviewed or obsolete accession. This version files it under the reviewed entry | 879 | 857 | 268 |
| `lysine_window` | the source window is centred on a lysine. The previous version turned it into a nearby arginine site | 132 | 71 | 0 |
| `truncated_window` | the window is cut short at a protein end. The previous version read it from the wrong end | 303 | 138 | 92 |
| `multi_arginine_peptide` | a peptide with several arginines that does not cover the reported position. The previous version guessed one | 39 | 5 | 0 |
| `peptide_start` | a peptide with several arginines, numbered differently from the entry. The previous version took the arginine nearest its middle. This version counts the reported position from the peptide's start in Maron 2021 Tables S2/S3 | 8 | 1 | 2 |
| `wrong_protein` | the previous version filed the position under another reviewed protein than the one whose sequence matched | 28 | 23 | 13 |
| `outdated_accession_or_symbol` | the source accession is merged or obsolete, or its gene symbol outdated. This version reaches the current entry | 103 | 0 | 30 |
| `evidence_mismatch` | the previous version kept the reported arginine although the source peptide or window does not match the entry there | 88 | 21 | 0 |
| `window_padding` | the window is padded where it runs off a protein end. The padding now gives the cut side | 5 | 0 | 1 |

A site can carry several codes when its records changed for different reasons.

## Proteins featured in the paper

| Protein | Sites in v1.0.0 | Sites in v1.1.0 | Removed | Added |
|---|---:|---:|---|---|
| FUS (P35637) | 30 | 34 | - | R422, R518, R522, R524 |
| EWSR1 (Q01844) | 35 | 36 | - | R518 |
| FMR1 (Q06787) | 13 | 13 | R562, R564 | R483, R558 |
| Histone H3.1 (P68431) | 10 | 10 | - | - |
| SmB/SNRPB (P14678) | 13 | 13 | - | - |
| SmD1/SNRPD1 (P62314) | 6 | 14 | - | R92, R104, R108, R110, R112, R114, R118, R119 |
| SmD3/SNRPD3 (P62318) | 6 | 7 | - | R125 |

## Decided in this version

- **Peptides with several arginines.** A peptide does not say which of its arginines is methylated, and the position the source reports may follow an older numbering. Where Maron et al. 2021 Tables S2 and S3 give the peptide's start in the same numbering, the reported position is counted from that start (reason `peptide_start`): 15 records (8 of them placed differently in v1.0.0), on ARID1B R608, ARID1B R640, ARID1B R648, CENPVL3 R17, TENT4A R145, WWC3 R493, WWC3 R500, WWC3 R507. Without such a start the record stays unplaced (reason `multi_arginine_peptide`). Where a start is known, counting from it agrees with all 2,949 placements the pipeline made from plain peptides by other means.
- **One UniProt reference, release 2026_03.** v1.0.0 used sequences staged on 2026-04-07 whose UniProt release was not recorded. Between them, 13 reviewed entries changed sequence; two carry sites: POLR2A (one substitution, outside every window; its sites keep their numbers) and RAB12, whose site R10 is R106 now.
- **Domain context for every protein.** The published analyses used InterPro intervals staged for only part of the union, so some proteins had no domain annotation. This version uses InterPro 108.0 for all reviewed entries: 453 proteins whose published sites were all "no domain annotation" now have a domain or repeat, and the domain context of 654 kept sites changes.
- **PRMT inhibitor calls** (Maron et al. 2021 Table S3). A peptide changes when its adjusted p value is at most 0.05. Each site, antibody and inhibitor gets a call over its peptides, and each site and inhibitor a call over its antibodies (`supporting/maron2021_inhibitors_v1.1.0.tsv`). The same combination gives the CARM1 knockout call over the two cell lines.

## Still open

- **Windows that match on one side only** (10 sites, flag `one_sided_window`; 10 of them in v1.0.0): the source window matches the entry exactly on one side of the reported arginine and not on the other, and no other reviewed entry holds it. Three of the ten windows occur in full in another UniProt isoform of the protein (GANAB R185, TCF7L2 R126, SNX12 R153; checked 2026-10-20), so they come from a splice variant. The pipeline leaves these records unplaced; the release keeps the sites, flagged and left out of every analysis, until the policy is decided: SQSTM1 R312, GANAB R185, GCC2 R303, RUSC2 R882, CTTNBP2 R451, TCF7L2 R126, SNX12 R153, DIS3 R144, DIS3 R147, DIS3 R150.
- **Paralog windows** (22 sites, flag `paralog_window`): a window matches a paralog exactly and the named entry within 2 residues. The records stay on the entry the source named, accepted in `config/remap_evidence_allowlist.tsv`.

## Tracing a change

- `supporting/site_changes_v1.0.0_to_v1.1.0.tsv`: one row per site of either version, with its change, reason codes, the records involved and, for a removed site, where those records are now (for an added site, where they were).
- `supporting/source_records_v1.1.0.tsv`: one row per source record, with what the study reported, the evidence sequence, where each version placed it, the independent evidence check and the reason code.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
