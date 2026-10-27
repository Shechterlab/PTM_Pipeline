# Changes from v1.4.0 to v1.4.1

| | v1.4.0 | v1.4.1 |
|---|---:|---:|
| Sites | 11,632 | 11,632 |
| Proteins | 4,439 | 4,439 |
| Sites kept pending | 0 | 0 |
| Source records | 67,516 | 67,516 |
| Source records placed | 38,325 | 38,325 |
| Studies and databases | 19 | 19 |

Of the 11,632 sites of v1.4.0 (pending ones included), 11,632 are kept and 0 removed; 0 sites were added. The records of v1.4.0 are read in the same order and all are placed as before. 

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## Protein-level tests within abundance deciles

Mass spectrometry finds abundant proteins most, and the methylated proteins are mostly abundant ones (median PaxDb percentile 74 for the 4,386 with a value), so a protein set rich in abundant proteins can look enriched among them whatever their methylation does. Each protein-level over-representation test is run again with its background split into PaxDb 6.1 abundance deciles (method `mantel_haenszel`, the background ending "each PaxDb abundance decile a stratum"; README.md, Statistics). The rows sit beside the tests they repeat, with the same sets and counts; no other value of the release changes. Both directions are tested, so a set can be depleted among proteins of similar abundance though enriched without the strata.

| Analysis | Sets tested | Enriched at q <= 0.05 | Of these, enriched within abundance deciles | Depleted within abundance deciles |
|---|---:|---:|---:|---:|
| PANTHER class groups | 11 | 5 | 3 | 4 |
| PANTHER protein classes | 25 | 8 | 6 | 10 |
| PANTHER subclasses | 127 | 31 | 14 | 20 |
| Functional classes | 13 | 10 | 9 | 3 |
| Functional subclasses | 27 | 22 | 20 | 3 |
| Reactome top-level pathways | 29 | 20 | 12 | 7 |
| Reactome pathways (10 to 500 proteins) | 1,705 | 597 | 205 | 63 |
| HPA locations | 48 | 9 | 7 | 6 |
| CD-CODE condensates, against other condensate proteins | 40 | 13 | 8 | 1 |
| CD-CODE condensates, against the proteome | 40 | 30 | 17 | 0 |
| DepMap dependency classes | 4 | 3 | 3 | 1 |
| GO RNA and DNA binding | 4 | 4 | 4 | 0 |
| Disease families (audited) | 8 | 8 | 2 | 1 |
| Disease families (additional) | 15 | 15 | 3 | 2 |

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (Site specificity in the CHANGES.md of v1.2.0): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.4.0_to_v1.4.1.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.4.1.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
