# Release v1.1.0: file specification

For the session working on `claude/epic-franklin-oni9ft`. Written 2026-10-20 from a review of `releases/v1.1.0/` at `4133e56`.

The website (`docs/` on `claude/awesome-thompson-7j9n92`, served at shechterlab.github.io/PTM_Pipeline) will read the release directly as its only source of site-level data, so that the table people download and the website can never disagree. For that, the release needs codes instead of prose in its columns, complete lists, and the numbers that are now packed into text. This document says what each file holds. **The union itself does not change.**

## Version

Rebuild as **v1.1.0** in place, unless v1.1.0 has already been deposited with a DOI (Zenodo or elsewhere). In that case build **v1.1.1**, and say in CHANGES.md that only the file format changed.

## Do not change

- The site set: 13,358 sites plus the 10 pending splice-junction sites. Also keep every placement, the evidence check, the reason codes and their meaning.
- `record_id` values. The other files cite them.
- The reference sequences (UniProt snapshot of 2026-04-07) and the annotations (InterPro 108.0, MobiDB export of 2026-03-30).

## Conventions for every TSV

- UTF-8, tab-separated, LF line endings, one header row, no quoting. No value contains a tab, a newline or `;`.
- A missing value is the empty string. Never write `NA`, `nan`, `None` or `-`.
- **Lists:** items are joined with `; ` (semicolon and one space), in the order the column states. Readers split on `;` and trim.
- **Booleans:** `true` or `false`.
- **Coded columns hold codes only.** Codes are lowercase snake_case, unless the code is a name such as `Rme1`, `ADMA` or `GSK591`. The wording of each code goes in `data_dictionary.tsv`.
- **Methylation:**
  - `Rme1` (monomethyl);
  - `Rme2` (dimethyl, type not stated);
  - `Rme2a` (asymmetric);
  - `Rme2s` (symmetric).

  Antibodies keep their names, always in this order: `MMA`, `ADMA`, `SDMA`.
- **Studies** are referred to by `study_id` (see `studies_v1.1.0.tsv`) in every file.
- **Positions** are 1-based integers on the reference sequence.

## Files

| File | Status | One row per |
|---|---|---|
| `human_methylarginine_sites_v1.1.0.tsv` | revised | site (including the 10 pending) |
| `human_methylarginine_proteins_v1.1.0.tsv` | **new** | protein with at least one row in the site table |
| `studies_v1.1.0.tsv` | **new** | study or database |
| `supporting/source_records_v1.1.0.tsv` | revised | source record |
| `supporting/maron2021_detection_v1.1.0.tsv` | **new** | site × Maron table × antibody |
| `supporting/maron2021_inhibitors_v1.1.0.tsv` | **new** | site × antibody × inhibitor |
| `supporting/site_changes_v1.0.0_to_v1.1.0.tsv` | revised | site of either version |
| `data_dictionary.tsv` | **new** | column of any TSV |
| `release.json` | **new** | the release (one object) |
| `human_methylarginine_sites_v1.1.0.xlsx` | revised | one sheet per TSV, plus README and dictionary sheets |
| `README.md`, `CHANGES.md`, `MD5SUMS` | revised | |

## 1. `human_methylarginine_sites_v1.1.0.tsv`

Sorted by `uniprot_accession`, then `position`. The "Change" column says what differs from the current file.

| Column | Type | Values and rule | Change |
|---|---|---|---|
| `site_id` | string | `{uniprot_accession}:R{position}` | |
| `status` | code | `site`, `pending` (the pending splice-junction rows) | **new** |
| `uniprot_accession` | string | reviewed canonical accession | |
| `uniprot_entry_name` | string | | |
| `gene` | string | primary gene name; empty when UniProt gives none (9 entries now). Never fill it with the entry name. | |
| `protein_name` | string | recommended name | |
| `position` | integer | | |
| `residue` | code | always `R` | |
| `window_15`, `window_31` | string | as now; `_` beyond a protein end | |
| `protein_length` | integer | | |
| `reference_sequence_md5` | string | as now (uppercase hex, UniProt format) | |
| `methyl_states_reported` | list<code> | `Rme1; Rme2`: degrees given by the source records (now ProMetheusDB only) | **codes** (was me1/me2) |
| `maron2021_s2_states` | list<code> | `Rme1; Rme2`: degree called in Maron 2021 Table S2 | **renamed** from `maron2021_ip_states`; codes |
| `maron2021_s3_antibodies` | list<code> | `MMA; ADMA; SDMA`, in that order: antibodies whose IP recovered the site in Table S3 | **renamed** from `maron2021_ip_antibodies` |
| `maron2021_GSK591_call` | code | `decreased`, `increased`, `mixed`, `no_change`; empty when not measured. Rule in section 6 | **renamed**; `no significant change` becomes `no_change` |
| `maron2021_MS023_call` | code | as for GSK591 | **renamed** |
| `motifs` | list<code> | in this order: `RG; RGG; GRG; GAR; RXR; CARM1_proline_rich; CARM1_hydrophobic_proline; DR; RD; ER; RE` | separator `; ` (was `, `) |
| `disorder` | code | `disordered`, `disorder_boundary`, `ordered`, `no_disorder_annotation` | **codes** |
| `disorder_edge_distance` | integer | as now (scripts/06); empty with `no_disorder_annotation` | |
| `domain` | code | `in_domain`, `boundary`, `inter_domain_linker`, `distal`, `no_domain_annotation` | **codes** |
| `nearest_domain_interpro` | string | InterPro accession of the nearest domain or repeat (the innermost when inside several) | **split** from `nearest_domain` |
| `nearest_domain_name` | string | its InterPro name | **split** |
| `nearest_domain_start`, `nearest_domain_end` | integer | its span on the reference sequence | **split** |
| `nearest_domain_class` | string | the class from `config/domain_class_ontology.json`, assigned as `scripts/07_summarize_domain_enrichment.py` `classify_residue` does; `No domain annotation` when there is none | **new** |
| `domain_side` | code | `inside`, `N_terminal_to_domain`, `C_terminal_to_domain`; empty when there is no domain (scripts/05 `nearest_domain_side`) | **new** |
| `domain_edge_distance` | integer | as now; 0 inside a domain | |
| `n_studies` | integer | | |
| `studies` | list<study_id> | in `studies_v1.1.0.tsv` order (most sites first) | **study_id codes** |
| `study_dois` | list<string> | same order as `studies` | |
| `source_tables` | list<code> | `maron2021_s5; prometheusdb_mmc2; ptm_master` (the tables its records come from) | **new** |
| `n_records` | integer | | |
| `record_ids` | list<integer> | ascending | separator `; ` (was `;`) |
| `reported_positions` | list<integer> | positions the records report, ascending; only when at least one differs from `position`, otherwise empty | **new** |
| `flags` | list<code> | in this order: `pending_splice_junction; paralog_window; near_match_only; shared_sequence; single_study` | separator `; ` |
| `same_window_accessions` | list<string> | **every** other reviewed entry whose sequence holds `window_31`, sorted. No "and N more". | **replaces** `same_window_in` |
| `same_window_genes` | list<string> | primary gene of each accession above, same order; empty item when there is none | **new** |
| `in_previous_version` | boolean | `true`/`false` | lowercase |
| `change_since_previous_version` | code | `kept`, `added` | |
| `change_reasons` | list<code> | reason codes (CHANGES.md) | separator `; ` |

**Removed:** `maron2021_inhibitor_detail`. Its numbers move to `supporting/maron2021_inhibitors_v1.1.0.tsv`.

Example, FUS R216 (selected columns):

```
site_id	status	gene	position	methyl_states_reported	maron2021_s2_states	maron2021_s3_antibodies	maron2021_GSK591_call	motifs	disorder	domain	nearest_domain_interpro	nearest_domain_name	nearest_domain_start	nearest_domain_end	domain_side	studies
P35637:R216	site	FUS	216	Rme1; Rme2	Rme1; Rme2	MMA; SDMA	no_change	RG; GRG; RXR	disordered	distal	IPR000504	RNA recognition motif domain	285	371	N_terminal_to_domain	Larsen2016; PSP; Maron2021; Fong2019; ProMetheusDB; Fedoriw2019; UniProt; Lim2020
```

## 2. `human_methylarginine_proteins_v1.1.0.tsv` (new)

Sorted by `uniprot_accession`. One row for each accession in the site table.

| Column | Type | Rule |
|---|---|---|
| `uniprot_accession` | string | |
| `uniprot_entry_name` | string | |
| `gene` | string | primary gene name; empty when there is none |
| `gene_synonyms` | list<string> | UniProt gene synonyms (the `synonym` rows of `inputs/reference/uniprot_accession_history_2026_03.tsv` for this accession), excluding `gene` |
| `reported_genes` | list<string> | gene symbols the source records placed on this protein used, excluding `gene` and `gene_synonyms`; sorted |
| `protein_name` | string | |
| `protein_length` | integer | |
| `reference_sequence_md5` | string | |
| `n_sites` | integer | rows with `status` = `site` |
| `n_pending` | integer | rows with `status` = `pending` |
| `studies` | list<study_id> | union over its sites, in studies-file order |

The website uses `gene_synonyms` and `reported_genes` as search aliases.

## 3. `studies_v1.1.0.tsv` (new)

Sorted by `n_sites`, descending. Every `study_id` used in any file appears here.

| Column | Type | Rule |
|---|---|---|
| `study_id` | code | see the list below |
| `label` | string | how people should see it |
| `source_table` | code | `maron2021_s5`, `prometheusdb_mmc2` or `ptm_master` |
| `citation` | string | short form, below |
| `doi` | string | |
| `licence` | string | e.g. `CC BY-NC-SA 3.0` for PSP; `see article` for journal supplements |
| `n_records`, `n_sites`, `n_sites_only_here` | integer | over `status` = `site` rows |
| `note` | string | PSP: `downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5`. Maron2021: `sites detected in Maron et al.'s own data (Table S5 reference "This study")` |

| study_id | label | citation | doi |
|---|---|---|---|
| `Larsen2016` | Larsen 2016 | Larsen SC et al., Sci Signal 9:rs9 (2016) | 10.1126/scisignal.aaf7329 |
| `PSP` | PhosphoSitePlus (PSP) | Hornbeck PV et al., Nucleic Acids Res 43:D512 (2015) | 10.1093/nar/gku1267 |
| `Maron2021` | Maron 2021 | Maron MI et al., iScience 24:102971 (2021) | 10.1016/j.isci.2021.102971 |
| `Guo2014` | Guo 2014 | Guo A et al., Mol Cell Proteomics 13:372 (2014) | 10.1074/mcp.O113.027870 |
| `Musiani2019` | Musiani 2019 | Musiani D et al., Sci Signal 12:eaat8388 (2019) | 10.1126/scisignal.aat8388 |
| `Fong2019` | Fong 2019 | Fong JY et al., Cancer Cell 36:194 (2019) | 10.1016/j.ccell.2019.07.003 |
| `ProMetheusDB` | ProMetheusDB | Massignani E et al., Mol Cell Proteomics 21:100243 (2022) | 10.1016/j.mcpro.2022.100243 |
| `Li2021` | Li 2021 | Li WJ et al., Nat Commun 12:1946 (2021) | 10.1038/s41467-021-21963-1 |
| `Fedoriw2019` | Fedoriw 2019 | Fedoriw A et al., Cancer Cell 36:100 (2019) | 10.1016/j.ccell.2019.05.014 |
| `UniProt` | UniProt | UniProt Consortium, Nucleic Acids Res 53:D609 (2025) | 10.1093/nar/gkae1010 |
| `Lim2020` | Lim 2020 | Lim Y et al., Proteome Sci 18:6 (2020) | 10.1186/s12953-020-00162-8 |
| `Wei2020` | Wei 2020 | Wei M et al., Front Cell Dev Biol 8:545934 (2020) | 10.3389/fcell.2020.545934 |
| `IEDB` | IEDB | Vita R et al., Nucleic Acids Res 47:D339 (2019) | 10.1093/nar/gky1006 |

The citations were checked against Crossref on 2026-10-20. The rows above are already in v1.1.0 site-count order. `PSP` as the code satisfies the PSP terms' requirement that bulk downloads name PSP as the source.

## 4. `supporting/source_records_v1.1.0.tsv` (revised)

Keep every column and its order. Change these columns:

| Column | Change |
|---|---|
| `source_table` | codes `maron2021_s5`, `prometheusdb_mmc2`, `ptm_master` |
| `study` | `study_id` code |
| `evidence_kind` | codes `window_13`, `window_15`, `window_31`, `peptide`, `marked_peptide`, `position_only` |
| `evidence_check` | codes: `agrees`, `agrees_no_site` (was "agrees: no site"), `one_of_several` ("one of several evidence sites"), `ambiguous` ("evidence ambiguous"), `dropped_evidence_has_site` ("dropped, evidence has a site"), `contradicts` ("contradicts evidence") |
| `change_since_previous_version` | codes `same`, `moved`, `newly_placed`, `no_longer_placed`, `not_placed` |
| `pending_site` | **new**: the `site_id` of the pending row this record supports; empty otherwise. The 10 pending rows' `record_ids` point at records with an empty `placed_site`; this column links them back. |

`placement` (words) and `placement_code` stay as they are.

## 5. `supporting/maron2021_detection_v1.1.0.tsv` (new)

One row per site × Maron table × antibody, for every row of the site table that Table S2 or S3 reaches.

| Column | Type | Rule |
|---|---|---|
| `site_id` | string | |
| `maron_table` | code | `S2` or `S3` |
| `antibody` | code | S3: `MMA`, `ADMA`, `SDMA`; S2: empty |
| `methyl_states` | list<code> | S2: `Rme1; Rme2` as called; S3: empty |
| `n_peptides` | integer | peptides of that table (and antibody) placed on the site |

## 6. `supporting/maron2021_inhibitors_v1.1.0.tsv` (new)

One row per site × antibody × inhibitor (Maron 2021 Table S3).

| Column | Type | Rule |
|---|---|---|
| `site_id` | string | |
| `antibody` | code | `MMA`, `ADMA`, `SDMA` |
| `inhibitor` | code | `GSK591` (PRMT5), `MS023` (type I PRMTs) |
| `n_peptides` | integer | peptides measured |
| `n_decreased` | integer | peptides with adjusted p ≤ 0.05 and log2 difference < 0 |
| `n_increased` | integer | peptides with adjusted p ≤ 0.05 and log2 difference > 0 |
| `log2_diff_median` | float | median `<inhibitor>_vs_DMSO_diff` over the peptides, 3 decimals |
| `log2_diff_strongest` | float | difference of the peptide with the lowest adjusted p, 3 decimals |
| `adj_p_min` | float | lowest `<inhibitor>_vs_DMSO_p.adj`, 3 significant digits |
| `call` | code | `decreased` when n_decreased > 0 and n_increased = 0; `increased` when the reverse; `mixed` when both > 0; `no_change` when both are 0 |

**Site-level call** (`maron2021_<inhibitor>_call` in the site table): combine the site's antibody rows for that inhibitor:

- `decreased` when at least one row is `decreased` and none is `increased` or `mixed`;
- `increased` likewise;
- `mixed` when rows disagree or any row is `mixed`;
- `no_change` when every row is `no_change`;
- empty when the site was not measured.

Significance is adjusted p **≤ 0.05**, as in scripts/16 (`--padj-cutoff`). The current release uses < 0.05; change it.

A peptide carrying several arginines gives each of them the same values. Say so in the dictionary.

## 7. `supporting/site_changes_v1.0.0_to_v1.1.0.tsv` (revised)

- `in_previous_version` and `in_this_version` become `true`/`false`.
- `reason_codes`, `issues`, `records_now_at_or_previously_at` and `record_ids` use the `; ` separator.
- `gene` is blank on 882 rows now, mostly removed sites on accessions that are not in the reference. Fill it from the records' `reported_gene` (most frequent) when the accession has no reference entry.

## 8. `data_dictionary.tsv` (new)

| Column | Rule |
|---|---|
| `file` | file name |
| `column` | column name |
| `type` | `string`, `integer`, `float`, `boolean`, `code`, `list<code>`, `list<integer>`, `list<string>` |
| `values` | for `code` and `list<code>`: every allowed code, joined with `; `; otherwise empty |
| `value_meanings` | for coded columns: `code = meaning` pairs joined with `; ` |
| `description` | one or two sentences, as in the README now |

Every column of every TSV appears exactly once, and nothing else does.

## 9. `release.json` (new)

```json
{
  "version": "v1.1.0",
  "previous_version": "v1.0.0",
  "release_date": "YYYY-MM-DD",
  "doi": "",
  "code_commit": "<40-character commit>",
  "build_time_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "counts": {"sites": 13358, "pending_sites": 10, "proteins": 4722, "proteins_with_pending": 4723,
             "source_records": 28960, "records_placed": 27983, "studies": 13,
             "kept": 12495, "removed": 1224, "added": 873},
  "reference": {"proteome": "UP000005640", "staged": "2026-04-07", "uniprot_release": "not recorded; before 2026_03",
                "accession_history": "2026_03"},
  "annotations": {"interpro": "108.0", "mobidb_export": "2026-03-30"},
  "phosphositeplus": {"download_date": "2021-05-22", "via": "Maron et al. 2021 Table S5",
                      "terms_archive": "inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf"},
  "licence": "CC BY-NC-SA 4.0",
  "files": [{"name": "human_methylarginine_sites_v1.1.0.tsv", "rows": 13368, "md5": "..."}],
  "inputs": [{"file": "sources/maron2021_tableS5_compiled_methylarginine_sites.xlsx", "md5": "..."}]
}
```

The counts must be computed, not typed in. `inputs` copies `inputs/MANIFEST.tsv`. The website reads this file for its version, counts and provenance.

## 10. README, CHANGES, workbook, checksums

- **Counts.** State them one way everywhere: "13,358 sites on 4,722 proteins, and 10 pending sites (`status` = `pending`), 13,368 rows". Now the README and CHANGES give 13,358 while the table has 13,368 rows.
- **PhosphoSitePlus.** State the download date: "PhosphoSitePlus data downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5; terms of use as archived on 2021-07-11 (inputs/licences/…)".
- **Columns section.** Generate it from `data_dictionary.tsv`, so the two can't drift.
- **Workbook.** One sheet per TSV, plus `README` and `data_dictionary` sheets. Each data sheet must equal its TSV.
- **MD5SUMS** lists every file in the folder.
- **v1.1.1 only:** if this becomes v1.1.1, CHANGES.md opens with "Format changes only; the sites, records and placements are those of v1.1.0" and lists the renamed and new columns.

## Checks the build must pass (stop on any failure)

1. `site_id` is unique and equals `{uniprot_accession}:R{position}`.
2. `residue` = `R` = `window_31`[16th character] = `window_15`[8th character]. Both windows match the reference sequence, and `reference_sequence_md5` matches it.
3. Rows with `status` = `site` plus rows with `status` = `pending` equal the row count. README, CHANGES and `release.json` give the same counts.
4. **Record links:**
   - every `record_ids` entry exists in `source_records`;
   - for `status` = `site`, every record's `placed_site` equals the `site_id`;
   - for `status` = `pending`, every record has an empty `placed_site` and `pending_site` equal to the `site_id`;
   - every placed record appears in exactly one site row.
5. Every study is in `studies_v1.1.0.tsv`. `n_studies` equals the length of `studies`; `study_dois` is aligned with `studies`; the counts in `studies_v1.1.0.tsv` recompute.
6. The accessions in `proteins` equal the site-table accessions; `n_sites` and `n_pending` recompute.
7. **Every coded column holds only its dictionary values.** The dictionary lists every column of every TSV and nothing else.
8. **No forbidden strings:** no `nan`, `None`, `NA`, `True`/`False` (capitalised), "and … more", or `, ` inside a list column.
9. **Maron tables:**
   - every `site_id` exists in the site table;
   - the site table's `maron2021_s3_antibodies` and `maron2021_s2_states` recompute from the detection table;
   - the two inhibitor calls recompute from the inhibitor table.
10. Every accession in `same_window_accessions` holds `window_31` in the reference FASTA. No other reviewed entry holds it.
11. MD5SUMS covers every file, and the workbook sheets equal the TSVs.
12. A second build from the same commit gives byte-identical TSVs.

## What the website will do with each file

| File | Website use |
|---|---|
| sites | every site, its windows, context, methylation evidence, flags (as badges and a filter) and "new in v1.1.0" |
| proteins | protein list, search aliases, lengths, sequence checksums (compared with live UniProt) |
| studies | study filter, "How to cite these sites", the PSP licence notice |
| source_records | the site card's reports: what each study reported, how it was placed, the independent check |
| maron2021_detection, maron2021_inhibitors | antibody and inhibitor panels and filters, per-antibody numbers on the site card |
| site_changes | an old link to a removed site (e.g. FMR1 R562) explains why it was removed |
| release.json | version, counts and the Dataset provenance section |

## After the release

The analyses the website shows are still those of v1.0.0 and must be re-run on this union, on `status` = `site` rows only:

- functional classes (scripts/03);
- domain, disorder and motif tests (05–08);
- other PTMs (09);
- clustering (10);
- condensates (12);
- cluster architecture (15, 25);
- the Maron summary (16);
- disease terms (34; needs R).

Before that, fix the re-run hazards listed in `DATA_LINEAGE_REVIEW.md`. The first: step 00 overwrites the staged InterPro intervals with empty files, and step 03 then stops.

The website export will refuse analysis tables whose site set differs from this release's `status` = `site` rows.
