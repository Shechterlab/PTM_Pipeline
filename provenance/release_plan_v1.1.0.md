# v1.1.0: release plan and file specification

For the release session (`claude/epic-franklin-oni9ft`) and the website session (`claude/awesome-thompson-7j9n92`). This is draft 2, of 2026-10-20. It replaces draft 1 (commit `f362340`), merges in the release session's proposed layout, and answers its two questions (section 9).

## Decisions

1. **One clean branch.** All of v1.1.0 is assembled on one new branch, `release-v1.1.0`, cut from `main`. The two working branches are archived as tags, not merged (section 1).
2. **The release is the single source of truth.** `releases/v1.1.0/` holds every table: sites, proteins, interactions, statistics. The website files in `docs/data/` are generated only from those tables. The only exception is what protein pages look up live in the browser (section 7).
3. **What goes in Git.**
   - Git holds the readable tables the build reads (uncompressed TSV, which Git packs about five-fold: the current 134 MB of `inputs/` is 28 MB packed).
   - Raw bulk downloads (`.gz`, `.zip`, `.xlsx` dumps, API responses) are attached to the `v1.1.0` GitHub release, and later to Zenodo with the paper's data DOI.
   - Each raw file is listed in `inputs/MANIFEST.tsv` with its URL, date, licence and MD5, and each tidy table names the raw file it came from.
4. **Interaction data is frozen, not downloaded again.**
   - The cache that built the current website is in the website session's container: `data/interactions_cache/`, 658 MB, 36,550 files. Its file dates are 2026-10-16, and 2026-10-19 for 6 files; the export that used it ran on 2026-10-19.
   - It also holds the pepR-MS supplement (`crosslinks/pepr_media-1.zip`), so nothing needs saving from bioRxiv by hand.
   - The website session turns it into tidy tables **first**, because that container is temporary.
5. **Keep every per-protein result the website shows:**
   - functional classes and subclasses;
   - arginine enrichment;
   - density clustering with its architecture classes;
   - disease families and terms;
   - CD-CODE condensates, HPA locations, other PTMs, the CARM1 knockout and interactions.
6. **Disease enrichment moves to Python.**
   - R runs once, to export the DOSE DisGeNET tables (`scripts/34_export_dgn_mapping.R`) and the org.Hs.eg.db Entrez-to-symbol table into `inputs/disease/`.
   - The per-term test (hypergeometric, Benjamini–Hochberg, terms of 10–300 genes) is then Python. It must reproduce the v1.0.0 `all_methyl_proteins_disease_enrichment.tsv` from the v1.0.0 union before it is used (section 6).
7. **File formats** are those of section 5. The union does not change: 13,358 sites and 10 pending, and every `record_id` as now. One site is renumbered by decision 8: RAB12 R10 becomes R106.
8. **One UniProt reference: release 2026_03** (released 2026-09-02), the reviewed human reference proteome UP000005640, in `inputs/reference/` on the website branch:
   - FASTA, entries, accession history and `uniprot_release.json`, with the release taken from UniProt's response header.
   - It replaces `uniprot_human_reviewed_2026-04-07.fasta`, whose release was never recorded and is older.
   - The accession history (already 2026_03), the interaction data and the layers are on 2026_03 already.
   - v1.1.0 moves to it. 13 sequences differ from the 2026-04-07 FASTA, and they carry 4 v1.1.0 sites:
     - POLR2A R426, R1603 and R1810 keep their positions;
     - RAB12 (Q6IQ22) gains 96 N-terminal residues, so its R10 becomes R106.
   - `inputs/reference/v1.1.0_sites_changed_in_2026_03.tsv` has the new values. Update `protein_length` and `reference_sequence_md5` for both proteins and every table that places the RAB12 site.
   - scripts/37 takes the same FASTA with `--reference-fasta`.
9. **Licences checked on 2026-10-20** (open decision 1 in section 9, now settled):
   - RBDmap: CC BY 4.0 (open access article, PMC5003815).
   - Trendel et al. 2025: CC BY-NC 4.0 (open access article).
   - CD-CODE: CC BY-NC 4.0 (cd-code.org).
   - iPTMnet: CC BY-NC-SA 4.0 (its licence page).
   - These and pepR-MS and DisGeNET are compatible with CC BY-NC-SA 4.0, so their tidy tables go in Git and their raw downloads in the public archive.
   - MobiDB states no licence on its own site; CC BY 4.0 is reported elsewhere. Publish its derived intervals with that note.
   - CPLM 4.0 is "All Rights Reserved" with no licence. Its download and tables are checksum-only, and the archive keeps them in a restricted file.

## 1. Branch plan

`main` (`36bc966`) is an ancestor of both working branches, so the clean branch can later fast-forward `main`.

1. **Freeze first** (website session, on its own branch). Write the tidy interaction tables and the staged layers (CPLM 4.0, HPA, CARM1 knockout) into `inputs/interactions/` and `inputs/layers/`, with MANIFEST rows. Record the MD5 of every cached file. If the PI wants the raw cache itself kept, the website session packs it as split archives and hands them over to attach to the GitHub release (it cannot upload release assets itself); otherwise those raw files are `checksum_only`.
2. **Cut `release-v1.1.0` from `main`.** Bring in the final files of both branches with `git checkout <branch> -- <paths>`, not a merge, so the branch history doesn't carry intermediate large files. Use one commit per area:
   1. code: `scripts/`, `config/`, `hpc/`, `envs/`, `tests/`;
   2. `inputs/`: tidy tables and MANIFEST;
   3. website code: `docs/index.html` and `docs/assets/`, without `docs/data/`.
3. **Remove from the tree:**
   - `PTM_results/` (v1.0.0 outputs, kept at tag `v1.0.0`);
   - `data/` and `results/`;
   - the April 2026 InterPro intervals, which only reproduce v1.0.0 and go to the v1.0.0 assets.

   Put what the change tracking needs (the v1.0.0 site table and source records) in `releases/v1.0.0/`.
4. **Build on the clean branch.** `bash hpc/build_release.sh v1.1.0` writes `releases/v1.1.0/`; then `scripts/36` and `37` write `docs/data/`. Nothing is copied from another branch's outputs.
5. **Preview.** Point GitHub Pages at `release-v1.1.0` `/docs` (Settings → Pages). Until then the live site stays on `claude/awesome-thompson-7j9n92`.
6. **Publish.**
   - Fast-forward `main` to `release-v1.1.0`, serve Pages from `main` `/docs`, and tag `v1.1.0`.
   - Attach the raw input parts to the GitHub release; deposit the release folder and the raw inputs on Zenodo.
7. **Archive.** Tag the working branches (`archive/awesome-thompson-2026-10`, `archive/epic-franklin-2026-10`), then delete the branches.

Each session pushes only its own branch. The PI creates `release-v1.1.0`, or says which session may push to it.

## 2. Layout

```
inputs/
  MANIFEST.tsv           one row per file: path, bytes, md5, source URL, version or download date, licence,
                         redistribution (git | release_asset | checksum_only), derived_from (for tidy tables)
  reference/             UniProt FASTA, entries, accession history
  sources/               Maron S2, S3, S5; ProMetheusDB; PTM master; licences/
  annotation/            InterPro 108.0 intervals, MobiDB intervals, CD-CODE tables
  interactions/          NEW: tidy tables frozen from the website cache (section 4)
  layers/                NEW: CPLM 4.0 lysine methylation, HPA locations, CARM1 knockout (tidy, as staged)
  disease/               NEW: DOSE DisGeNET term-gene table, Entrez-symbol table (exported from R once)
config/                  ontologies, allow-list
scripts/                 pipeline; build_release.py; 36 and 37 read releases/ only
releases/
  v1.0.0/                published union: sites and source records (for change tracking)
  v1.1.0/                section 5
docs/                    website; docs/data/ generated from releases/v1.1.0/
```

## 3. Inputs

**Rules for every input file:**
- Named `<source>_<content>_<version or YYYY-MM-DD>.<ext>`.
- Listed in `MANIFEST.tsv`, with checksums verified before every build (as `hpc/build_release.sh` does now).
- **A tidy table** (filtered or reshaped from a raw download) names its raw file in `derived_from`, and is written by a script in `scripts/` that runs from the raw file.
- **A raw file not in Git** has `redistribution` = `release_asset` (attached to the GitHub release) or `checksum_only` (when its licence doesn't allow redistribution; the tidy table may still be in Git, see section 9).

**What exists and stays:**
- the reference sequences;
- the three source tables and Maron S2/S3;
- InterPro 108.0 and MobiDB intervals;
- the CD-CODE tables;
- the accession history;
- the PTM master.

## 4. Inputs to freeze from the website cache

Written by a new `scripts/37a_stage_interactions.py`, the parsing half of today's `scripts/37`. It runs once on the cache; download dates are the cache file dates.

| Tidy table (`inputs/interactions/`) | From (cache) | Licence |
|---|---|---|
| `string_v12.5_physical_links_ge400.tsv` | `9606.protein.physical.links.v12.5.txt.gz` | CC BY 4.0 |
| `complexportal_complexes_<date>.tsv` (curated and hu.MAP 3.0 predicted, one row per complex and member) | `9606.tsv`, `9606_predicted.tsv` | CC0 |
| `quickgo_rna_dna_binding_<date>.tsv` | `quickgo_rna.tsv`, `quickgo_dna.tsv` | CC BY 4.0 |
| `uniprot_reader_writer_domains_<date>.tsv` (PRMT, Tudor, WD40, PHD) | `uniprot_human_*.tsv` | CC BY 4.0 |
| `uniprot_dna_binding_regions_<date>.tsv` | `crosslinks/uniprot_dna_bind.tsv` | CC BY 4.0 |
| `pdbe_interface_residues_<date>.tsv` (accession, residue, partner kind, PDB entries) | `pdbe_interfaces/` | EMBL-EBI terms |
| `sifts_observed_segments_<date>.tsv` | `uniprot_segments_observed.tsv.gz`, `pdb_chain_uniprot.tsv.gz` | EMBL-EBI terms |
| `pdb_site_contacts_<date>.tsv` (residues within 4 Å of each interface methylarginine) | `contacts_seqres/` | CC0 |
| `pdb_interface_environments_<date>.tsv` (solvent accessibility, Levy class, base or backbone H-bonds) | `environments/` | CC0 |
| `pdb_methylarginine_residues_<date>.tsv` (2MR, DA2, NMM) | `modelled_methyl_r/` | CC0 |
| `protvar_interface_residues_2024.05.28.tsv` (accession, residue, partner, pDockQ) | `2024.05.28_interface_summary_5A.tsv.gz` | CC BY 4.0 |
| `predicted_model_contacts_<date>.tsv` (measured salt bridges and H-bonds; ProtVar served models and AlphaFold DB complexes) | `protvar_measure/`, `afdb_measure/`, `afdb_complex/`, `af2complexes_interfaces_v1.tsv.gz` | CC BY 4.0 |
| `intact_features_<release>.tsv` (human participants, sequence-checked) | `intact/` | CC BY 4.0 |
| `pepr_ms_rna_crosslinks_v2.tsv` | `crosslinks/pepr_media-1.zip` | CC BY-NC |
| `rbdmap_fragments_castello2016.tsv` | `crosslinks/rbdmap_castello2016_table_s1.xlsx` | ORA terms; check |
| `trendel2025_dna_crosslinks.tsv` | `crosslinks/trendel2025_mmc1.xlsx`, `mmc2.xlsx` | Elsevier supplement; check |
| `cdcode_confidence_api2.3_<date>.tsv` | `cdcode/` | not stated; check |

**Layers** (`inputs/layers/`), from the staged files in the same container:
- `cplm4_lys_methyl_validated.tsv` (and `_rejected`; CPLM licence not stated);
- `hpa_subcellular_locations_v25.1.tsv` (CC BY 4.0);
- `carm1_ko_shishkova2017_peptides.tsv` and `_sites.tsv` (CC BY 4.0).

**Disease** (`inputs/disease/`, exported from R once and recorded with the package versions):
- `dose_<version>_disgenet_term_gene.tsv`;
- `orghsegdb_<version>_entrez_symbol.tsv`.

## 5. Release files (`releases/v1.1.0/`)

### Conventions for every TSV

- UTF-8, tab-separated, LF, one header row, no quoting; no value contains a tab, newline or `;`.
- Missing is the empty string, never `NA`, `nan`, `None` or `-`.
- **Lists:** `; `-separated, in the order the column states.
- **Booleans:** `true`/`false`.
- **Coded columns** hold codes only (lowercase snake_case, or a name such as `Rme1`, `ADMA`, `GSK591`); meanings go in `data_dictionary.tsv`.
- **Methylation:** `Rme1`, `Rme2` (type not stated), `Rme2a`, `Rme2s`.
- **Antibodies:** `MMA`, `ADMA`, `SDMA`, in that order.
- **Studies** by `study_id`. **Positions** 1-based on the reference sequence.

### Files

| File | Status | One row per |
|---|---|---|
| `human_methylarginine_sites_v1.1.0.tsv` | revised | site (including the 10 pending) |
| `human_methylarginine_proteins_v1.1.0.tsv` | **new** | protein |
| `studies_v1.1.0.tsv` | **new** | study or database |
| `site_interactions_v1.1.0.tsv` | **new** | site × piece of interaction evidence |
| `protein_partners_v1.1.0.tsv` | **new** | protein × STRING physical partner (score ≥ 0.7) |
| `protein_complexes_v1.1.0.tsv` | **new** | protein × complex |
| `statistics_v1.1.0.tsv` | **new** | statistical test the website or the paper reports |
| `supporting/source_records_v1.1.0.tsv` | revised | source record |
| `supporting/maron2021_detection_v1.1.0.tsv` | **new** | site × Maron table × antibody |
| `supporting/maron2021_inhibitors_v1.1.0.tsv` | **new** | site × antibody × inhibitor |
| `supporting/carm1_ko_v1.1.0.tsv` | **new** | site × cell line |
| `supporting/disease_terms_v1.1.0.tsv` | **new** | enriched DisGeNET term |
| `supporting/residue_evidence_v1.1.0.tsv` | **new** | protein residue with interface, crosslink or predicted-interface evidence (all residues, for protein maps) |
| `supporting/other_ptm_sites_v1.1.0.tsv` | **new** | other modification on a protein of the union |
| `supporting/site_changes_v1.0.0_to_v1.1.0.tsv` | revised | site of either version |
| `data_dictionary.tsv` | **new** | column of any TSV |
| `release.json` | **new** | the release |
| `human_methylarginine_sites_v1.1.0.xlsx` | revised | sheets: README, sites, proteins, studies, data_dictionary |
| `README.md`, `CHANGES.md`, `MD5SUMS` | revised | |

### 5.1 Sites

Sorted by `uniprot_accession`, then `position`. "Change" compares with the current file.

| Column | Type | Values and rule | Change |
|---|---|---|---|
| `site_id` | string | `{uniprot_accession}:R{position}` | |
| `status` | code | `site`, `pending` | **new** |
| `uniprot_accession`, `uniprot_entry_name` | string | | |
| `gene` | string | primary gene name; empty when UniProt gives none (9 entries); never the entry name | |
| `protein_name` | string | | |
| `position` | integer | | |
| `residue` | code | always `R` | |
| `window_15`, `window_31` | string | `_` beyond a protein end | |
| `protein_length` | integer | | |
| `reference_sequence_md5` | string | uppercase hex, UniProt format | |
| `methyl_states_reported` | list<code> | `Rme1; Rme2` from the source records | **codes** (was me1/me2) |
| `maron2021_s2_states` | list<code> | `Rme1; Rme2` (Table S2) | **renamed** from `maron2021_ip_states` |
| `maron2021_s3_antibodies` | list<code> | `MMA; ADMA; SDMA` (Table S3) | **renamed** from `maron2021_ip_antibodies` |
| `maron2021_GSK591_call`, `maron2021_MS023_call` | code | `decreased`, `increased`, `mixed`, `no_change`; empty when not measured (rule in 5.6) | **renamed**; new rule |
| `carm1_ko_call` | code | as above, over cell lines (5.7) | **new** |
| `motifs` | list<code> | order: `RG; RGG; GRG; GAR; RXR; CARM1_proline_rich; CARM1_hydrophobic_proline; DR; RD; ER; RE` | separator |
| `disorder` | code | `disordered`, `disorder_boundary`, `ordered`, `no_disorder_annotation` | **codes** |
| `disorder_edge_distance` | integer | scripts/06; empty with no annotation | |
| `domain` | code | `in_domain`, `boundary`, `inter_domain_linker`, `distal`, `no_domain_annotation` | **codes** |
| `nearest_domain_interpro`, `nearest_domain_name` | string | | **split** from `nearest_domain` |
| `nearest_domain_start`, `nearest_domain_end` | integer | | **split** |
| `nearest_domain_class` | string | `config/domain_class_ontology.json` as `scripts/07` `classify_residue` assigns it; `No domain annotation` when none | **new** |
| `domain_side` | code | `inside`, `N_terminal_to_domain`, `C_terminal_to_domain`; empty when none (scripts/05) | **new** |
| `domain_edge_distance` | integer | 0 inside | |
| `interaction_status` | code | `experimental`, `predicted_only`, `none` (rule in 5.4) | **new** |
| `interaction_evidence` | list<code> | evidence types at the site, in the order of 5.4 | **new** |
| `n_studies` | integer | | |
| `studies` | list<study_id> | studies-file order | **codes** |
| `study_dois` | list<string> | same order | |
| `source_tables` | list<code> | `maron2021_s5; prometheusdb_mmc2; ptm_master` | **new** |
| `n_records` | integer | | |
| `record_ids` | list<integer> | ascending | separator |
| `reported_positions` | list<integer> | only when one differs from `position` | **new** |
| `flags` | list<code> | order: `pending_splice_junction; paralog_window; near_match_only; shared_sequence; single_study` | separator |
| `same_window_accessions` | list<string> | every other reviewed entry holding `window_31`, sorted; never truncated | **replaces** `same_window_in` |
| `same_window_genes` | list<string> | primary gene of each, same order | **new** |
| `in_previous_version` | boolean | | lowercase |
| `change_since_previous_version` | code | `kept`, `added` | |
| `change_reasons` | list<code> | CHANGES.md codes | separator |

Removed: `maron2021_inhibitor_detail`; its numbers move to `maron2021_inhibitors`.

Example (selected columns): `P35637:R216  site  FUS  216  Rme1; Rme2  Rme1; Rme2  MMA; SDMA  no_change  RG; GRG; RXR  disordered  distal  IPR000504  RNA recognition motif domain  285  371  N_terminal_to_domain  Larsen2016; PSP; Maron2021; Fong2019; ProMetheusDB; Fedoriw2019; UniProt; Lim2020`.

### 5.2 Proteins

Sorted by accession; one row for each accession in the site table.

| Group | Columns |
|---|---|
| Identity | `uniprot_accession`, `uniprot_entry_name`, `gene`, `gene_synonyms` (UniProt synonyms from the accession history), `reported_genes` (other symbols the records used), `protein_name`, `protein_length`, `reference_sequence_md5` |
| Sites | `n_sites`, `n_pending`, `studies`, `n_sites_experimental`, `n_sites_predicted_only` |
| Function | `functional_classes`, `functional_subclasses` (lab ontology, scripts/03), `go_rna_binding`, `go_dna_binding` (code: `experimental`, `other_evidence` or empty; QuickGO GO:0003723 and GO:0003677 with child terms; experimental = EXP, IDA, IPI, IMP, IGI, IEP and their high-throughput forms) |
| Arginine enrichment | `arginine_count`, `arg_enrichment_log2` (shrinkage estimate, scripts/08), `arg_enrichment_q` (empty when not tested) |
| Clustering (scripts/10, 15, 25) | `nearest_spacing_min`, `nearest_spacing_median`, `max_sites_3aa`, `max_sites_5aa`, `max_sites_10aa`, `max_sites_20aa`, `max_sites_50aa`, `best_window_aa`, `best_window_start`, `best_window_end`, `best_window_sites`, `best_window_log2_obs_null`, `best_window_q`, `dense_windows` (list of `start-end`, q ≤ 0.05), `architecture_class` |
| Context | `condensates` (CD-CODE ids), `hpa_locations` (`location (reliability)`), `disease_families`, `als_ftd_audited` (boolean) |
| Interactions | `n_string_partners` (≥ 0.7), `n_partners_methylated`, `n_complexes` |
| Other PTMs | `n_phospho`, `n_lys_acyl`, `n_ub_sumo`, `n_lys_methyl` (positions in `other_ptm_sites`) |

### 5.3 Studies

Columns:
- `study_id`, `label`, `source_table`, `citation`, `doi`, `licence`;
- `n_records`, `n_sites`, `n_sites_only_here`;
- `note`, with the PSP note: "downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5".

The ids, labels and Crossref-checked citations, in site-count order:

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

### 5.4 Site interactions

One row per site and piece of evidence. Columns:

| Column | Rule |
|---|---|
| `site_id` | |
| `evidence_type` | `rna_crosslink`, `rbdmap_fragment`, `dna_crosslink`, `pdb_protein`, `pdb_rna`, `pdb_dna`, `pdb_self`, `predicted_interface`, `predicted_contact`, `intact_mutation`, `intact_binding_region`, `intact_methylated` (this order is the order of `interaction_evidence` in 5.1) |
| `localisation` | `at_residue`, `candidate_residue`, `within_3`, `in_fragment`, `in_peptide`, `in_region` |
| `partner_kind` | `protein`, `rna`, `dna`, `self` |
| `partner_accession`, `partner_gene` | for protein partners |
| `source` | `pepR-MS`, `Kramer2014`, `XRNAX`, `RBS-ID`, `pRBS-ID`, `RBDmap`, `Trendel2025`, `PDBe-KB`, `ProtVar`, `AFDB`, `IntAct` |
| `structures` | list of PDB ids or model ids |
| `score` | pDockQ or ipTM for models |
| `distance_a` | for measured contacts |
| `detail` | short text (e.g. the IntAct effect, the bond type) |

Rule for `interaction_status`, as on the website now:
- **experimental:** any row of `rna_crosslink` (at the residue or as a candidate residue), `rbdmap_fragment`, `dna_crosslink` (at or candidate), or `pdb_*`;
- **predicted_only:** otherwise, any `predicted_*` row;
- **none:** neither.

IntAct rows never set the status.

### 5.5 Partners, complexes and residue evidence

- **`protein_partners`:** `uniprot_accession`, `partner_accession`, `partner_gene`, `string_score`, `partner_in_union`, `partner_classes` (`prmt; tudor; wd40; phd`).
- **`protein_complexes`:** `uniprot_accession`, `complex_id`, `complex_name`, `curation` (`curated`, `predicted`), `n_members`, `n_members_methylated`.
- **`residue_evidence`:** every residue of a union protein with evidence, for the protein maps.
  - Columns: `uniprot_accession`, `position`, `evidence_type` (codes of 5.4), `n_partners`, `best_score`, `sources`.
  - It replaces today's `methylarg_pred_iface.json` and the residue lists inside `methylarg_interactions.json`.

### 5.6 Maron 2021 tables

- **`maron2021_detection`:** `site_id`, `maron_table` (`S2`, `S3`), `antibody` (S3 only), `methyl_states` (S2 only), `n_peptides`.
- **`maron2021_inhibitors`:** `site_id`, `antibody`, `inhibitor` (`GSK591`, `MS023`), `n_peptides`, `n_decreased`, `n_increased`, `log2_diff_median`, `log2_diff_strongest`, `adj_p_min`, `call`.

**Rule:**
- A peptide is significant at adjusted p **≤ 0.05** (scripts/16 `--padj-cutoff`).
- Per row: `decreased` when it has decreased and no increased peptides, `increased` the reverse, `mixed` both, `no_change` neither.
- Site call (per inhibitor, over its antibody rows): `decreased` when at least one row is decreased and none is increased or mixed; `increased` likewise; `mixed` when rows disagree or any is mixed; `no_change` when all are.
- Peptides with several arginines give each the same values.

### 5.7 CARM1 knockout, disease terms, other PTMs

- **`carm1_ko`:** `site_id`, `cell_line` (`MCF7`, `MDA-MB-231`), `n_peptides`, `log2_ko_over_wt`, `p`, `call`. The rule is that of scripts/41: |log2| ≥ 1 and p < 0.05. Combine over cell lines as in 5.6.
- **`disease_terms`:** `term_id`, `term_name`, `family`, `n_genes_term`, `n_methylated`, `p`, `p_adjust`, `genes` (symbols). Terms of 10–300 genes with p_adjust ≤ 0.05; families by the audited rules of scripts/34.
- **`other_ptm_sites`:** `uniprot_accession`, `position`, `residue`, `ptm_type` (`phospho`, `lys_acyl`, `ub_sumo`, `lys_methyl`), `source` (`ptm_master`, `cplm4`), `evidence`.
  - The PTM master rows on union proteins are included, not only half of them as the website shows now (see `DATA_LINEAGE_REVIEW.md`).
  - CPLM rows only if their licence allows (section 9).

### 5.8 Statistics

One row per test that the website or the paper reports. The Overview charts and the Interactions findings read only this table.

| Column | Rule |
|---|---|
| `analysis` | code, e.g. `functional_class`, `functional_subclass`, `disease_family`, `motif`, `disorder`, `domain_context`, `domain_class`, `condensate_conf3`, `hpa_location`, `rna_crosslink`, `dna_crosslink_candidate`, `pdb_interface_rna`, `pdb_interface_protein`, `predicted_interface`, `intact_binding_region`, `go_rna_binding`, `partners_methylated`, `complex` |
| `category` | the bar or row (e.g. `RNA-binding / RNP`) |
| `stratum` | e.g. `all`, `ordered_only`, `in_domain`, `no_rg` |
| `a`, `b`, `c`, `d` | the 2×2 counts (target with, target without, background with, background without) |
| `odds_ratio`, `ci_low`, `ci_high`, `p`, `q` | |
| `method` | `fisher`, `mantel_haenszel`, `hypergeometric` |
| `background` | what was compared (e.g. `other arginines of the same proteins`) |
| `n_strata` | for Mantel–Haenszel |

### 5.9 Source records, site changes, dictionary, release.json

- **`source_records`:**
  - columns as now, with codes for `source_table`, `study` (study_id), `evidence_kind` (`window_13`, `window_15`, `window_31`, `peptide`, `marked_peptide`, `position_only`), `evidence_check` (`agrees`, `agrees_no_site`, `one_of_several`, `ambiguous`, `dropped_evidence_has_site`, `contradicts`) and `change_since_previous_version` (`same`, `moved`, `newly_placed`, `no_longer_placed`, `not_placed`);
  - a new `pending_site` column.
- **`site_changes`:** booleans `true`/`false`, lists `; `. Fill the 882 blank `gene` values from the records' `reported_gene`.
- **`data_dictionary.tsv`:** `file`, `column`, `type`, `values` (allowed codes), `value_meanings`, `description`. It lists every column of every TSV, and nothing else.
- **`release.json`:**
  - version, previous version, date, DOI, code commit, build time;
  - computed counts;
  - reference proteome (UP000005640, UniProt release 2026_03; `inputs/reference/uniprot_release.json`);
  - annotation versions (InterPro 108.0, MobiDB 2026-03-30, STRING 12.5, CD-CODE 2.3, HPA 25.1, IntAct release);
  - PhosphoSitePlus (download 2021-05-22, via Maron Table S5, terms archive);
  - licence;
  - files with rows and MD5;
  - inputs copied from MANIFEST.
- **README and CHANGES:**
  - State the counts one way everywhere: "13,358 sites on 4,722 proteins, and 10 pending sites, 13,368 rows".
  - Generate the column sections from the dictionary.

## 6. Analyses to compute on v1.1.0

Run on `status` = `site` rows only. **Before** a step runs on v1.1.0, it must reproduce its v1.0.0 output from the v1.0.0 union and inputs; a step that differs stops the build. This catches changes from the re-run itself, as distinct from changes from the union.

| Result | Script | Inputs | Check against v1.0.0 |
|---|---|---|---|
| Functional classes and enrichment | 03 | PTM master, InterPro, `config/functional_ontology.json` | class tables identical |
| Domain and disorder context, enrichment | 05, 06, 07 | InterPro 108.0, MobiDB | identical with the April intervals; then 108.0 (658 kept sites change, as CHANGES.md says) |
| Motifs and arginine enrichment | 08 | FASTA | identical |
| Other PTMs | 09 | PTM master | identical |
| Clustering and architecture | 10, 15, 25 | union, FASTA | identical with the same seed and permutations |
| Condensates | 12, 37b | CD-CODE tables and API 2.3 table | identical |
| Maron summary | 16 | Maron S2/S3 | identical |
| Disease terms and families | new Python step, 34 | `inputs/disease/` | same terms, counts and p.adjust (to 1e-9) as the R run |
| Interactions | 37b | `inputs/interactions/` | the v1.0.0 site export gives today's `methylarg_interactions.json` content |
| CARM1, HPA, CPLM layers | 39, 41, 42 merge | `inputs/layers/` | identical |

`DATA_LINEAGE_REVIEW.md` lists re-run hazards that must be fixed first. Among them: step 00 overwrites the InterPro intervals with empty files; one plot reads the published tables; several steps fetch from moving URLs. Every step reads `inputs/`, never the network.

## 7. Website generation

- **`scripts/36`** reads only `releases/v1.1.0/`. It writes `docs/data/methylarg_sites.json` and `methylarg_evidence.json` as compact encodings of the release tables: vocabularies, bit masks, arrays.
- **`scripts/37b`** (the computing half of today's 37) reads `inputs/interactions/` and the release sites, and writes the interaction tables of section 5. `docs/data/methylarg_interactions.json` is then encoded from those tables.
- **A round-trip check** decodes each JSON and compares it with the tables. The build fails if any site, protein, statistic or label differs, or if a JSON value has no table behind it.
- **What stays live** in the browser, as now, documented on the Documentation page:
  - UniProt features;
  - AlphaFold DB models and complexes;
  - PDBe and RCSB structure lists and files;
  - STRING partners below 0.7 and links among partners;
  - ENCODE, ChIP-Atlas, Open Targets, Ensembl Compara, NCBI titles.

## 8. Checks

Draft 1's twelve checks stay:
1. unique `site_id`;
2. residue, windows and MD5 match the FASTA;
3. counts agree everywhere;
4. record links;
5. studies;
6. proteins;
7. codes only from the dictionary;
8. no forbidden strings;
9. the Maron calls recompute;
10. complete same-window lists;
11. MD5SUMS and workbook;
12. byte-identical rebuild.

New checks:

13. `interaction_status` and `interaction_evidence` recompute from `site_interactions`.
14. Every protein in `protein_partners`, `protein_complexes`, `residue_evidence` and `other_ptm_sites` is in the proteins table.
15. Every `analysis` code in `statistics` is in the dictionary, and its counts recompute from the release tables.
16. The website round-trip of section 7.
17. Every file under `inputs/` matches MANIFEST, and every tidy table regenerates byte-identically from its `derived_from` file when that file is available.

## 9. Answers to the release session's two questions, and open decisions

**Which per-protein statistics to keep:** all four:
- functional classes;
- arginine enrichment;
- density clustering with the architecture classes;
- disease families.

The website uses each in filters, protein pages and Overview charts. All but disease are Python and cheap. Disease needs the one-time R export above, then Python.

**70 MB more raw downloads in Git?** No.
- Raw downloads (mostly already compressed, so Git can't shrink them) go to the GitHub release assets and Zenodo.
- Git keeps the tidy uncompressed tables, which Git packs well. The current `inputs/` is 134 MB checked out but 28 MB packed. The new tidy interaction tables should add a few tens of MB checked out, and far less packed.

**Open, for the PI:**
1. ~~Redistribution of sources whose licence is unclear or restrictive.~~ Settled (decision 9).

**Archive.** The website branch's `zenodo/` folder describes the offline and Zenodo archive of the raw downloads:
- a public part (every source whose licence allows it);
- a restricted part (CPLM);
- a Git bundle of every branch;
- checksums, deposit metadata and an upload script.

The raw downloads go there rather than into GitHub release assets.
2. **Pending splice-junction sites on the website:** shown flagged, or hidden.
3. **Inhibitor rule:** confirm adjusted p ≤ 0.05 and the per-antibody combination of 5.6.
4. **Version:** rebuild as v1.1.0 in place, unless v1.1.0 already has a DOI; then v1.1.1.
5. **Branch:** who creates `release-v1.1.0`, and which session pushes to it.
