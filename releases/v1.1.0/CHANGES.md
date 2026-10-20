# Changes from v1.0.0 to v1.1.0

v1.1.0 uses the sources and methods of DeAngelo et al. (2026), mapped to UniProt release 2026_03, with corrections to how some source records were placed on UniProt sequences. v1.0.0 is the published union (git tag v1.0.0, `releases/v1.0.0/`).

| | v1.0.0 | v1.1.0 |
|---|---:|---:|
| Sites | 13,719 | 13,365, and 6 pending |
| Proteins | 5,090 | 4,723 |
| Source records placed | 28,134 of 28,964 | 27,996 of 28,964 |

Of the 13,719 published sites, 12,495 are kept (sites of v1.1.0 too, 4 of them pending) and 1,224 were removed; 876 sites were added. A site moves when its records now point at another arginine: the old one is removed and the new one added. A kept site can still have gained or lost records: 2,384 did, and at 31 none of the previous records remain, other records now naming the same arginine (supporting/site_changes, column records_change). One kept site has a new number only: RAB12 R10 is R106 in UniProt 2026_03, which adds 96 residues at its N-terminus.

## What was corrected

Each record whose placement changed carries one reason code.

| Reason code | What changed | Records | Sites removed | Sites added |
|---|---|---:|---:|---:|
| `peptide_recentred` | the previous version moved the peptide to the arginine nearest its middle. This version keeps the reported arginine inside the peptide | 1,773 | 117 | 482 |
| `unreviewed_accession` | the previous version filed reviewed-entry numbering under the source's unreviewed or obsolete accession. This version files it under the reviewed entry | 879 | 857 | 268 |
| `lysine_window` | the source window is centred on a lysine. The previous version turned it into a nearby arginine site | 132 | 71 | 0 |
| `truncated_window` | the window is cut short at a protein end. The previous version read it as if centred, taking the arginine at or nearest the middle of the short window. This version counts from its intact end | 303 | 138 | 92 |
| `multi_arginine_peptide` | a peptide with several arginines that does not cover the reported position. The previous version guessed one | 35 | 3 | 0 |
| `peptide_start` | a peptide with several arginines, numbered differently from the entry. The previous version took the arginine nearest its middle. This version counts the reported position from the peptide's start in Maron 2021 Tables S2/S3 | 8 | 1 | 2 |
| `wrong_protein` | the previous version filed the position under another reviewed protein than the one whose sequence matched | 28 | 23 | 13 |
| `outdated_accession_or_symbol` | the source accession is merged or obsolete, or its gene symbol outdated. This version reaches the current entry | 88 | 0 | 23 |
| `entry_not_in_reference` | no reviewed entry has the source accession or gene symbol (an entry deleted from Swiss-Prot, a readthrough or an unreviewed entry). The evidence occurs once in the reviewed proteome, and this version places it there. For the deleted duplicate genes CBSL, OVOS2 and SIK1B that is the paralog holding the same sequence (CBS, OVOS1, SIK1) | 15 | 0 | 7 |
| `evidence_mismatch` | the previous version kept the reported arginine although the source peptide or window does not match the entry there | 81 | 23 | 2 |
| `window_padding` | the window is padded where it runs off a protein end. The padding now gives the cut side | 5 | 0 | 1 |

A site can carry several codes when its records changed for different reasons.

## Corrections to the source tables

- **Windows missing NAN** (14 records, `source_correction` `nan_restored`). Maron et al. 2021 Table S5 gives these Larsen 2016 windows 28 residues long, where every other Larsen window has 31, and no Larsen window there holds the tripeptide NAN: the letters were read as a missing value when the table was made. The parser restores NAN where inserting it at exactly one place gives the window centred on the reported position of the reported entry; the window as the table gives it is kept in `evidence_as_reported`. v1.0.0 placed 14 of these records; this version places 14.
- **Records without an accession** (4 records, `source_correction` `accession_na`). Table S5 gives them a gene symbol and a window, and the accession NA; v1.0.0 left them out. They are placed from their gene symbol and window (2 placed) and numbered after the records of v1.0.0, so that every other record keeps its record_id.

## Proteins featured in the paper

| Protein | Sites in v1.0.0 | Sites in v1.1.0 | Removed | Added |
|---|---:|---:|---|---|
| FUS (P35637) | 30 | 34 | - | R422, R518, R522, R524 |
| EWSR1 (Q01844) | 35 | 36 | - | R518 |
| FMR1 (Q06787) | 13 | 13, and 1 pending | R562, R564 | R483, R558, R606 (pending) |
| Histone H3.1 (P68431) | 10 | 10 | - | - |
| SmB/SNRPB (P14678) | 13 | 13 | - | - |
| SmD1/SNRPD1 (P62314) | 6 | 14 | - | R92, R104, R108, R110, R112, R114, R118, R119 |
| SmD3/SNRPD3 (P62318) | 6 | 7 | - | R125 |

## Decided in this version

- **Peptides with several arginines.** A peptide does not say which of its arginines is methylated, and the position the source reports may follow an older numbering. Where Maron et al. 2021 Tables S2 and S3 give the peptide's start in the same numbering, the reported position is counted from that start (reason `peptide_start`): 15 records (8 of them placed differently in v1.0.0), on ARID1B R608, ARID1B R640, ARID1B R648, CENPVL3 R17, TENT4A R145, WWC3 R493, WWC3 R500, WWC3 R507. Without such a start the record stays unplaced (reason `multi_arginine_peptide`). Where a start is known, counting from it agrees with all 2,949 placements the pipeline made from plain peptides by other means.
- **One UniProt reference, release 2026_03.** v1.0.0 used sequences staged on 2026-04-07 whose UniProt release was not recorded. Between them, 13 reviewed entries changed sequence; two carry sites: POLR2A (one substitution, outside every window; its sites keep their numbers) and RAB12, whose site R10 is R106 now.
- **Domain context for every protein.** The published analyses used InterPro intervals staged for only part of the union, so some proteins had no domain annotation. This version uses InterPro 108.0 for all reviewed entries: 454 proteins whose published sites were all "no domain annotation" now have a domain or repeat, and the domain context of 656 kept sites changes.
- **Protein classes from PANTHER.** The paper's functional classes are keyword classes that overlap and leave much of the proteome out. This version adds the PANTHER 19.0 protein classes (scripts/03b): each protein's top-level class, grouped for the overview (config/panther_class_groups.json) so that the groups and the proteins PANTHER leaves unclassified partition the reviewed proteome, with tests of the groups, the 24 top-level classes and their subclasses (statistics: panther_group, panther_class, panther_subclass). The paper's classes are kept.
- **Disease families.** The ALS / motor neuron / FTD rule now also takes TDP-43 proteinopathy, Pick disease (not Niemann-Pick) and anterior horn cell disease, which the published rule missed (5 DisGeNET terms, among them the enriched terms TDP-43 Proteinopathies, Pick Disease of the Brain and Degeneration of anterior horn cells); the other 7 audited families are unchanged. The DisGeNET terms none of them takes now go to 15 additional, rule-based families (scripts/disease_families.py ADDITIONAL_FAMILY_RULES: mostly clinical signs, by organ system), tested and BH-adjusted apart (statistics: disease_family_additional). They are not audited.
- **Methylation degree at each site** (`methyl_state`, 9,967 of 13,365 sites). v1.0.0 gave a degree only from ProMetheusDB. This version also counts Maron et al. 2021 Table S2, UniProt's modified-residue annotations with experimental evidence (inferred ones are listed apart) and the two studies that looked for one degree only: Larsen 2016 (monomethylation) and Wei 2020 (ADMA antibody, dimethylation). Mass spectrometry tells mono- from dimethylation, but asymmetric from symmetric mostly by the antibody used, so `methyl_state` gives Rme1, Rme2 or both, and the columns it is drawn from keep the detail: 8,230 Rme1, 1,049 Rme1 + Rme2, 688 Rme2.
- **PRMT inhibitor calls** (Maron et al. 2021 Table S3). A peptide changes when its adjusted p value is at most 0.05. Each site, antibody and inhibitor gets a call over its peptides, and each site and inhibitor a call over its antibodies (`supporting/maron2021_inhibitors_v1.1.0.tsv`). The same combination gives the CARM1 knockout call over the two cell lines.

## Still open

- **Windows that match on one side only** (6 sites, flag `one_sided_window`; 4 of them in v1.0.0): the source window matches the entry exactly on one side of the arginine and not on the other, and no other reviewed entry holds it. Three of these windows occur in full in another UniProt isoform of the protein (GANAB R185, TCF7L2 R126, SNX12 R153; checked 2026-10-20), so they come from a splice variant. Some are reported at another position, in the numbering of an isoform, and match one complete flank (15 residues and the arginine) here: FMR1 R606 (reported at R564), MTCL1 R1269 (reported at R303). The pipeline leaves these records unplaced; the release keeps the sites, flagged and left out of every analysis, until the policy is decided: FMR1 R606, SQSTM1 R312, GANAB R185, TCF7L2 R126, SNX12 R153, MTCL1 R1269.
- **Paralog windows** (22 sites, flag `paralog_window`): a window or peptide matches another reviewed entry exactly and the named entry within 2 residues. The other entry is a paralog, except at RBM39 R109, where the 6-residue peptide SALRGK reads SAIRGK in RBM39 (I and L have the same mass) and also occurs in CCDC13. The records stay on the entry the source named, accepted in `config/remap_evidence_allowlist.tsv`.
- **Arginine picked by the reported position** (859 sites, flag `arginine_from_position`): every record of the site has a peptide with several arginines, or a window or peptide found more than once in the entry, so the evidence alone does not single out this arginine; the positions the sources report do.

## Tracing a change

- `supporting/site_changes_v1.0.0_to_v1.1.0.tsv`: one row per site of either version, with its change, reason codes, the records involved and, for a removed site, where those records are now (for an added site, where they were; for a kept site, the records that moved onto or off it, and both).
- `supporting/source_records_v1.1.0.tsv`: one row per source record, with what the study reported, the evidence sequence, where each version placed it, the independent evidence check and the reason code.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
