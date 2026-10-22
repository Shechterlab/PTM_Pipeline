# Changes from v1.1.0 to v1.2.0

| | v1.1.0 | v1.2.0 |
|---|---:|---:|
| Sites | 13,365 | 11,601 |
| Proteins | 4,723 | 4,432 |
| Sites kept pending | 6 | 0 |
| Source records | 28,964 | 66,660 |
| Source records placed | 27,996 | 37,518 |
| Studies and databases | 13 | 18 |

Of the 13,371 sites of v1.1.0 (pending ones included), 11,201 are kept and 2,170 removed; 400 sites were added. The records of v1.1.0 are read in the same order and 20,228 of them are placed differently (supporting/source_records, columns change_since_previous_version and reason_code). The records read for the first time (each study's own table) follow them, from record_id 28,965.

## Site specificity

A record adds a site only where its source localises the methyl group to one arginine: a localisation score over the study's threshold (MaxQuant localisation probability, ptmRS, Ascore or ModScore, ModLS), a peptide whose mass leaves one arrangement (no other residue the search could methylate, or every such residue a dimethylated arginine), a UniProt annotation with experimental evidence from the literature, or a PhosphoSitePlus site, trusted as curated (`localised_by` in `supporting/source_records`). Its peptide, or without one its window, must also name exactly one arginine in the reviewed proteome, in a protein the record names.

Every study is read from its own tables (`inputs/sources/<study>_methylarginine_records.tsv`, staged by `scripts/00_stage_study_sites.py`), so the rows of the compiled tables give way to them: 18,991 records of Maron et al. 2021 Table S5 and the PTM master table are not placed (`superseded_by_study_table`), the study's own records standing in for them. 8,487 records carry no localisation at all (`unresolved_not_localised`): 2,105 Fedoriw 2019, 1,486 Li 2021, 1,244 Larsen 2016, 1,236 ProMetheusDB, 588 Lu 2025, 506 Maron 2021, 324 Musiani 2019, 322 Fong 2019, 235 Geoghegan 2015, 175 Sylvestersen 2014, 169 Guo 2014, 59 Lim 2020, 38 IEDB. All records stay in `supporting/source_records`.

The 2,170 sites of v1.1.0 that no localised record supports are removed. By the sources that reported them in v1.1.0:

| Reported in v1.1.0 by | Sites removed |
|---|---:|
| Larsen 2016 (with or without PhosphoSitePlus) | 1,086 |
| PhosphoSitePlus only | 373 |
| Maron 2021 (with or without PhosphoSitePlus) | 240 |
| Li 2021 (with or without PhosphoSitePlus) | 126 |
| two or more studies | 64 |
| ProMetheusDB (with or without PhosphoSitePlus) | 64 |
| Fedoriw 2019 (with or without PhosphoSitePlus) | 53 |
| Musiani 2019 (with or without PhosphoSitePlus) | 46 |
| Guo 2014 (with or without PhosphoSitePlus) | 36 |
| Lim 2020 (with or without PhosphoSitePlus) | 34 |
| Fong 2019 (with or without PhosphoSitePlus) | 26 |
| IEDB (with or without PhosphoSitePlus) | 16 |
| Wei 2020 (with or without PhosphoSitePlus) | 3 |
| UniProt (with or without PhosphoSitePlus) | 2 |
| Larsen 2016 and Maron 2021 (with or without PhosphoSitePlus) | 1 |

47 of them have a localised record whose peptide or window names several arginines in the reviewed proteome (paralogs, repeats) or one in a protein the record does not name.

The others were reported through the compiled tables by studies whose own tables do not localise them: the study's table gives the site below its localisation threshold, gives no score for a peptide that leaves several residues the search could methylate, or does not list it. Their records stay in `supporting/source_records` (`superseded_by_study_table`, or `unresolved_not_localised` for ProMetheusDB).

## Lower-confidence sites

2,406 arginines that a study reports as methylated, or that v1.0.0 or v1.1.0 had as sites, are not sites of this release: each of their records fails at least one site test. They are kept as lower-confidence sites (`human_methylarginine_lower_confidence_sites_v1.2.0.tsv`), with the results of the record that passes most tests, and are never in this release's counts and tests. The analyses were run a second time with the 1,844 whose own evidence covers the arginine and leaves the methyl group on an arginine (`statistics_with_lower_confidence_v1.2.0.tsv`, `analyses_with_lower_confidence/`). The website shows them, hatched, when the reader switches Lower-confidence sites on, and its counts and charts then come from that second run.

Every record has passed its study's identification filter (false discovery rate or score) before it is read. It is then put to five tests at the arginine its study reports; a site needs a record that passes all five.

| Test | Passes when | Lower-confidence sites whose best record fails it |
|---|---|---:|
| `evidence_at_site` | the record's own peptide or window occurs in this entry and covers this arginine, where the study puts the methyl group (for a peptide with several arginines: one of them) | 213 |
| `named_entry` | the arginine is in a protein the record names: its accession, the current entry of a merged accession, or an entry of its gene (not another protein that alone holds the sequence). Under the paralog rule, the record names one of the copies that hold the sequence | 10 |
| `unique_in_proteome` | the peptide or window names one site in the reviewed human proteome (UniProt 2026_03, UP000005640): no paralog, repeat or other protein holds it (for a peptide with several arginines: the peptide occurs once). Result paralogs: it names one site in each of several entries, one of them the protein the record names (the paralog rule, which counts as passing: the site goes on the copy most records of that sequence name) | 259 |
| `methyl_on_arginine` | the methyl group can only be on an arginine: the record is localised to one, or its peptide holds no lysine (a lysine could carry the methyl group instead), untested when the study gives a window but no peptide | 351 |
| `localised` | the study localises the methyl group to this arginine (localised_by: a localisation score over its threshold, a peptide whose mass leaves one arrangement, or a curated annotation) | 1,957 |

What became of every site of the earlier releases (`supporting/published_sites`):

| Release | Sites | `site` | `lower_confidence` | `not_in_reference` |
|---|---:|---:|---:|---:|
| v1.0.0 | 13,719 | 11,651 | 2,051 | 17 |
| v1.1.0 | 13,371 | 11,201 | 2,170 | 0 |

`site`, `lower_confidence`, `not_in_reference`: a site of this release; a lower-confidence site of this release (human_methylarginine_lower_confidence_sites); the accession is not a reviewed entry of the reference (an unreviewed or retired entry), or the position is not an arginine of it, and its records now name no arginine of the reference.

## Lower-confidence evidence

10,151 records add no site but say where a methyl group could be: they are not localised, or localised but not to one site of a protein they name. They are listed with their candidate sites in `supporting/lower_confidence_evidence`, never counted as sites or used in a test, and the website shows them only on request (protein map, row Lower-confidence evidence).

| Why | Records |
|---|---:|
| `several_residues`: the methyl group is on one of the candidate arginines of the peptide (or a lysine), not localised | 5,323 |
| `below_threshold`: the study's localisation score is below its threshold: the candidate is the residue the study names | 2,478 |
| `not_in_reference`: localised, but its peptide or window is not in the reviewed proteome (an isoform or another sequence) | 1,278 |
| `no_localising_evidence`: the study names the residue, but no spectrum or score it gives localises it | 638 |
| `several_sites`: localised, but its peptide or window occurs at several sites (paralogs, repeats): one of the candidates | 245 |
| `other_protein`: localised, but its peptide or window names one site, in a protein the record does not name | 141 |
| `search_disagrees`: the study's confident spectra disagree on which candidate carries the methyl group | 48 |

## Studies read from their own tables

The studies compiled in Maron et al. 2021 Table S5, and UniProt and IEDB, read from their own tables. Records: those the staging keeps (localised, and within the study's false discovery rate); sites: those they support now, beside the sites the study supported in v1.1.0 through the compiled tables.

| Study | Records | Placed | Sites | Sites in v1.1.0 |
|---|---:|---:|---:|---:|
| [Larsen 2016](https://doi.org/10.1126/scisignal.aaf7329) | 8,882 | 7,297 | 7,260 | 8,507 |
| [Maron 2021](https://doi.org/10.1016/j.isci.2021.102971) | 4,279 | 3,749 | 1,482 | 1,821 |
| [Guo 2014](https://doi.org/10.1074/mcp.O113.027870) | 3,258 | 2,999 | 1,569 | 1,444 |
| [Musiani 2019](https://doi.org/10.1126/scisignal.aat8388) | 2,511 | 2,093 | 1,299 | 1,395 |
| [Fong 2019](https://doi.org/10.1016/j.ccell.2019.07.003) | 2,543 | 2,185 | 1,270 | 1,297 |
| [Li 2021](https://doi.org/10.1038/s41467-021-21963-1) | 1,821 | 329 | 319 | 1,055 |
| [Fedoriw 2019](https://doi.org/10.1016/j.ccell.2019.05.014) | 2,659 | 536 | 409 | 925 |
| [UniProt](https://doi.org/10.1093/nar/gkae1010) | 254 | 254 | 212 | 739 |
| [Lim 2020](https://doi.org/10.1186/s12953-020-00162-8) | 1,350 | 1,254 | 561 | 603 |
| [Wei 2020](https://doi.org/10.3389/fcell.2020.545934) | 289 | 283 | 283 | 279 |
| [IEDB](https://doi.org/10.1093/nar/gky1006) | 81 | 43 | 43 | 84 |

`STUDIES.md` gives, for each, the tables read, the thresholds, and records traced from the source table to the reviewed proteome.

## Studies added (5)

Each study's supplementary table is staged like the others: one record per reported peptide and arginine, kept when it passes the study's false discovery rate and the methyl group is localised to that arginine, and placed by the rule above.

| Study | Sample | Records | Placed | Sites | Sites new to the union | Proteins new to the union |
|---|---|---:|---:|---:|---:|---:|
| [Hartel 2019](https://doi.org/10.1074/mcp.RA119.001625) | HEK293T; LN229 | 2,806 | 2,684 | 1,447 | 108 | 24 |
| [Geoghegan 2015](https://doi.org/10.1038/ncomms7758) | Jurkat; PBMCs | 4,536 | 4,069 | 2,233 | 168 | 45 |
| [Sylvestersen 2014](https://doi.org/10.1074/mcp.O113.032748) | HEK293T | 1,281 | 993 | 889 | 20 | 10 |
| [Lu 2025](https://doi.org/10.1038/s44319-025-00599-y) | HeLa | 743 | 155 | 155 | 10 | 0 |
| [Ma 2023](https://doi.org/10.1021/jasms.3c00154) | breast cancer PDX | 403 | 403 | 403 | 44 | 16 |

400 sites are new to the union, on 310 proteins (120 of them new).

### ADMA and SDMA from neutral losses

Hartel et al. 2019 tell asymmetric from symmetric dimethylarginine by the neutral loss in the fragment spectra: dimethylamine (45.058 Da) for ADMA, methylamine (31.042 Da) for SDMA. `dimethyl_type_neutral_loss` gives the call at 179 sites: 140 ADMA, 29 SDMA, 10 ADMA; SDMA. Maron et al. 2021 Table S3 says which antibodies recovered a site (`maron2021_s3_antibodies`). At the 90 sites with both a neutral-loss call and an ADMA or SDMA immunoprecipitation, the two agree at 24, overlap at 54 and differ at 12. Both columns are kept as reported.

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## Methylation degree

`methyl_state` gives a degree at 10,352 of 11,601 sites: 7,950 Rme1, 1,433 Rme1 + Rme2, 969 Rme2. The studies read from their own tables count through their records (`methyl_states_reported`).

## Protein-set tests

- Every test of protein or gene sets (PANTHER and functional classes, CD-CODE condensates, Human Protein Atlas locations, disease families, Gene Ontology binding, complexes) is an over-representation analysis with clusterProfiler::enricher, and the DisGeNET terms with DOSE::enrichDGN (scripts/enrich_sets.R; statistics: method enricher, with fold enrichment): p is the one-sided hypergeometric probability, fold enrichment (k/n)/(M/N), q Benjamini-Hochberg over the sets of the analysis (README, Statistics). The disease families are tested against the reviewed human genes with a DisGeNET term.
- The Human Disease Ontology is added, with the gene sets of DOSE enrichDO (`supporting/disease_ontology`), and every DisGeNET term tested is listed (`supporting/disease_terms_tested`).
- The analyses run on the sites of this version.

## Corrections to the analyses

- Functional classes: short gene tokens of `config/functional_ontology.json` matched inside other words (`atr` in matrix and tetratricopeptide, `actin` in interacting, `lamin` in filamin, `nop` in aminopeptidase); they now match whole words or word starts. Every protein is labelled from the same text (UniProt and InterPro); the names in the source tables, which list mostly methylated proteins, are no longer added to theirs. The classes called enriched do not change.
- Condensates: the background holds reviewed entries only, secondary accessions counted as their current entry.
- Disease enrichment: a protein counts under its UniProt primary gene (the source tables' first symbol could be a neighbouring gene, e.g. IGF2 for IGF2-AS), and the Disease Ontology table names genes by symbol.
- Density windows: a window of W residues spans W residues (it spanned W+1), each protein's null is seeded from its accession, and the null takes 10,000 permutations (1,000 before).
- Motifs: the proteome background leaves out the methylated arginines. Per-protein arginine enrichment: q over the proteins tested.
- Disorder: a protein end is not an edge of a disordered region; MobiDB entries numbered on another sequence than the reference were replaced from the MobiDB API or dropped (`inputs/annotation/mobidb_entries_patched_2026-10-21.tsv`).
- Human Protein Atlas: an accession several HPA genes map to takes the best reliability and the main location of any.
- Protein names keep brackets that belong to the name (DNA (cytosine-5)-methyltransferase 1); gene synonyms are split at the semicolons between genes.
- Interaction evidence: predicted contacts and interface residues whose residues disagree with the reference sequence (a chain of another protein, another numbering) and ProtVar pairs with a non-human partner are left out; PDB contacts in fusion chains name the protein of the segment that holds the residue; CD-CODE members listed under isoforms count once, with their best confidence.

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (above): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.1.0_to_v1.2.0.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.2.0.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
