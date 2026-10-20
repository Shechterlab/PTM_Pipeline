# Data lineage review (2026-10-19)

Where every input of the methylarginine analysis and its website comes from, what is frozen, what is not, and what can silently change results. Written for review before any change to `main`. Code reviewed: `claude/awesome-thompson-7j9n92` at `8e96f2f` (the branch that serves the website), plus the frozen inputs of `inputs-april-2026-update` (`b0bdc5f`). Checksums of the frozen inputs are in `inputs/MANIFEST.tsv`.

## Summary

- **The frozen inputs are the published run's inputs.** The three source tables are byte-identical to the copies committed on 2026-04-07 (`cd56dd0`). The UniProt FASTA in `inputs/` is the snapshot the published union was mapped to: the pre-fix remapping code run on it reproduces the published remapping on 28,960 of 28,960 rows (13,719 sites on 5,090 accessions), and all 12,846 saved sequence windows match it. Its UniProt release was never recorded; it predates 2026_03 (13 sequences differ, among them RAB12, 244 residues then, 340 now).
- **On that snapshot the corrected union is stable.** The 2026-10-18 code gives 13,324 sites on 4,706 proteins, the same as on UniProt 2026_03; the current code (with the accession history, issue 7) gives 13,353 on 4,720: 103 more rows placed, none moved or lost. Both pass the evidence check (22 accepted issue-9 rows) and the residue audit.
- **The PTM master has no build code anywhere.** It is a human iPTMnet-format export, joined to protein annotations, added as a finished file. Its release and download date are unknown. It holds no PhosphoSitePlus data, so the other-PTM layer is far from exhaustive.
- **The website shows about half of the master's validated other-PTM sites**, because of an export filter (33,203 of 69,669 dropped; 17,289 of 37,145 phospho sites on methylarginine proteins shown). The paper's cross-PTM analyses use the full set.
- **A rerun as the scripts stand would fail or mix data.** The first pipeline step overwrites the staged InterPro intervals with empty files and step 03 then stops; one extended plot reads the published tables instead of the new ones; several inputs are fetched at run time from moving targets. Details under Hazards.
- **The website's data is not reproducible from the repository.** Nothing that `scripts/37` downloaded (STRING, Complex Portal, PDBe-KB, crosslink sets, ProtVar, IntAct, CD-CODE API) or that scripts 39–42 staged (CPLM 4.0, CARM1 knockout, Human Protein Atlas) is frozen in any branch, and the deployed export was never rebuilt with the provenance fields.

## 1. Frozen inputs (`inputs/`)

| File (former name) | What it is | Status |
|---|---|---|
| `sources/ptm_master_iptmnet_human.tsv` (`human_ptm_master.tsv`) | Cross-PTM master (section 3) | identical to `cd56dd0`; origin undocumented |
| `sources/maron2021_tableS5_compiled_methylarginine_sites.xlsx` (`maron2021_TableS5_Compiled_Methylarginine_Data.xlsx`) | Maron et al. 2021 Table S5 | identical to `cd56dd0`; same cells as the publisher's mmc6.xlsx (re-saved, so the bytes differ) |
| `sources/maron2021_tableS2_methylarginine_ip_sites.xlsx` (new) | Maron et al. 2021 Table S2 (PTMScan IP sites, me1/me2) | the publisher's mmc3.xlsx; `scripts/16` on it reproduces `PTM_results/maron_state_treatment/` |
| `sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx` (new) | Maron et al. 2021 Table S3 (ADMA/SDMA/MMA IP after GSK591 or MS023) | the publisher's mmc4.xlsx; same check |
| `sources/massignani2022_prometheusdb_mmc2_methylsites.xlsx` (`promteheusdb_supplement_mmc2.xlsx`) | ProMetheusDB mmc2 (Massignani et al. 2022) | identical to `cd56dd0` and byte-identical to the publisher's mmc2.xlsx |
| `reference/uniprot_human_reviewed_2026-04-07.fasta`, `_entries.tsv` (`uniprot_human_reviewed_canonical.fasta`, `_metadata.tsv`) | Reference sequences, staged 2026-04-07 | proven to be the published run's snapshot; release unknown |
| `reference/uniprot_accession_history_2026_03.tsv` (`uniprot_accession_history.tsv`) | Secondary accessions and gene synonyms, UniProt 2026_03 | built 2026-10-19 against the reference FASTA |
| `annotation/mobidb_disorder_intervals_2026-03-30.tsv` (`mobidb_human_reviewed_disorder_intervals.tsv`) | Disorder intervals from the MobiDB export of 2026-03-30 | covers all reviewed proteins; the source JSON is not frozen |
| `annotation/interpro_domain_intervals_2026-04-07.tsv` (`interpro_human_reviewed_domain_like_intervals.tsv`) | Domain intervals from InterPro `protein2ipr.dat.gz` (current release on 2026-04-07) | **covers only the published union's 5,090 accessions**; 81 proteins of the corrected union (111 sites), 94 with the current code (126 sites), are missing; release number unknown |
| `annotation/interpro_all_entry_intervals_2026-04-07.tsv` (`interpro_human_reviewed_all_entry_intervals.tsv`) | All InterPro entries, same accessions | read by no script |
| `annotation/cdcode_protein_membership_v2.2.tsv`, `cdcode_proteins_2026-03-18.csv`, `cdcode_condensates_2026-03-18.csv` (`protein2cdcode_v2.2.tsv`, `proteins_202603181653.csv`, `condensates_202603181647.csv`) | CD-CODE v2.2 membership and the 2026-03-18 export tables | version and date only in the file names |

`inputs/MANIFEST.tsv` gives each file's checksum, former name, source URL, version, citation, licence and how it was checked. Only `hpc/build_release.sh` reads `inputs/`; `hpc/ptm_env.sh` still points at `PTM_DATA_ROOT/{base,external,context,condensates}` with the old file names.

## 2. Inputs that are not frozen

| Input | Used by | Affects | How it is obtained now |
|---|---|---|---|
| UniProt accession history | `scripts/02` (issue 7) | union (+103 rows) | downloaded by step 00 from UniProt `current_release` `sec_ac.txt` when missing, so it can be newer than the FASTA (2026_03 here) |
| InterPro `protein2ipr.dat.gz` (16 GB), `entry.list` | step 00 rebuild of the intervals | domain context, functional classes | your local copy of 2026-04-07; otherwise downloaded from `current_release` |
| MobiDB `mobidb_human_2026-03-30.json` (611 MB) | `00_stage_mobidb_disorder.py`, which no driver calls | only to rebuild the frozen intervals | your local copy |
| Maron Table S2/S3 | `scripts/16` (methyl state, IP, inhibitor response) | website annotations | nothing stages them; the rerun copies the published `maron_state_treatment/` tables |
| CPLM 4.0 methylation, CARM1 knockout table, Human Protein Atlas | scripts 39–42 | website layers and `scripts/40` | downloaded from moving URLs; checksums recorded in the export, files not kept |
| Everything behind `methylarg_interactions.json` | `scripts/37` | website interactions | cached in `data/interactions_cache`, which Git ignores and no branch holds |
| GO, DisGeNET, org.Hs.eg.db | scripts 13, 14, 34 | GO and disease enrichment | R packages; `envs/ptm_pipeline_hpc.yml` pins minor versions only (DOSE 4.4.*, org.Hs.eg.db 3.22.*; GO.db and HDO.db unpinned); no run records `sessionInfo()` |

## 3. The PTM master

No code builds it on any branch or in any commit (including the `scripts.zip` once committed). Early READMEs call it an "existing iPTMnet/UniProt-based master". Its columns follow iPTMnet's `ptm.txt` plus a protein-annotation join and derived columns.

| Source code | Rows | Database |
|---|---|---|
| hprd | 64,133 | HPRD (phospho) |
| unip | 42,559 | UniProtKB |
| glygen | 25,012 | GlyGen (glyco) |
| iedb | 9,092 | IEDB |
| rlim, rlim+ | 6,984, 429 | RLIMS-P text mining |
| sign | 5,508 | SIGNOR |
| npro | 5,107 | neXtProt (inferred) |
| pro | 2,917 | Protein Ontology |
| sno | 2,362 | dbSNO |
| intact | 823 | IntAct |

- 164,926 rows, all human: phosphorylation 119,602; glycosylation 29,478; methylation 7,395; acetylation 5,124; S-nitrosylation 2,362; ubiquitination and SUMOylation 698.
- Arginine methylation: 847 rows (762 UniProt, 85 IEDB); lysine methylation: 245 rows. Most "methylation" rows are IEDB methylated cysteines (6,146), which get no PTM group.
- 51,083 rows (31%, mostly HPRD) carry isoform accessions (`-1` in 32,412); on all of them the annotation columns are blank.
- No PhosphoSitePlus, Phospho.ELM, dbPTM or CPLM rows. iPTMnet itself integrates PhosphoSitePlus, so either its bulk download omits it or it was filtered out; which is not recorded.
- Newest PubMed ID about late 2023, so the export is probably from 2024 or later; the iPTMnet release is unknown.

**Why the website looks short on other PTMs.** Two separate reasons:

1. Coverage: the master is a modest public set (no PhosphoSitePlus; ubiquitin, acetylation and lysine methylation are one to two orders of magnitude below proteomics-scale resources). The website now adds CPLM 4.0 for lysine methylation: of the 3,259 methyl-K entries it shows, 179 come from the master.
2. An export bug. `scripts/09_compare_ptm_contexts.py:66` keeps the first row of each duplicated site, which in file order is often HPRD's `-1` row; `scripts/36_export_site_explorer_data.py:451` then keeps only rows whose accession equals the canonical one. Of 69,669 sequence-validated other-PTM sites, 33,203 are dropped (33,096 phospho). 14,071 of them also have a canonical-accession row in the master and are lost only through the order of rows; the other ~19,100 are isoform-numbered rows that passed a single-residue check, and need a window check before they can be shown. The paper's cross-PTM statistics (scripts 09 and 11) use the full validated set and are not affected.

**What would document or rebuild it.** A builder that downloads iPTMnet's `ptm.txt` and protein file, keeps human rows, maps PTM groups as the frozen table does, maps `-1` accessions to canonical with a sequence check, and records release, date and checksum; tested by reproducing the frozen table. Whether to add PhosphoSitePlus (licence) or dbPTM is a scientific choice that changes the cross-PTM results.

## 4. Hazards (ranked)

1. **Step 00 empties the staged InterPro intervals.** `hpc/slurm_run_integration.sh:81` runs `scripts/00_stage_context_sources.py --skip-interpro` on `PTM_CONTEXT_ROOT`; with no InterPro rows the script still writes `interpro_human_reviewed_*_intervals.tsv` (the domain file becomes a single newline) and rewrites `context_source_summary.tsv`. Step 03 (`slurm_run_integration.sh:132-136`) then raises `EmptyDataError` before InterPro is restaged (lines 183-195). Reproduced on a copy. The rerun script and `ISSUES_2026-10-19.md` (issue 17) both say the staged inputs are reused unchanged; they are not. Never point `PTM_CONTEXT_ROOT` at `inputs/`.
2. **InterPro intervals exist only for the run's accessions** (restaged with `--accessions` of the new union, `slurm_run_integration.sh:185,195`). Domain context and the InterPro text behind the functional classes exist only for methylated proteins; the frozen file lacks the proteins that enter the corrected union. Rebuild once for all reviewed proteins from the same `protein2ipr.dat.gz` and freeze that file.
3. **Step 24 plots the published tables.** `scripts/run_extended_analyses.sh:84-85` passes only `--outdir`; `scripts/24_plot_domain_edge_length_null.py:146-151` defaults to `PTM_results/domain_context/...` relative to the working directory, so a rerun from the repository plots the published tables into the new results root. Scripts 33–35 and the R scripts 13–14 also default to the repository's `PTM_results`; the drivers override all but script 24.
4. **Run-time downloads from moving targets:** UniProt REST and `sec_ac.txt` (`current_release`), InterPro `current_release` and API, Human Protein Atlas, CPLM, CD-CODE `releases/latest`, Complex Portal `current`, IntAct `current`, QuickGO, PDBe-KB, SIFTS, RCSB, AlphaFold DB. Only STRING (12.5) and ProtVar (2024.05.28) are pinned. If step 00 finds no FASTA in `PTM_CONTEXT_ROOT` it downloads today's UniProt even when `PTM_CANONICAL_FASTA` points elsewhere.
5. **Silent skips that change results:** no accession history (issue 7 rows stay unplaced), no MobiDB intervals (disorder-matched nulls and disorder context dropped), no CD-CODE tables (non-human and synthetic memberships kept), no `Rscript` (GO and disease enrichment skipped), missing staged website inputs.
6. **No run record in the main pipeline.** Only the rerun writes `RUN_INFO.txt`, and it hashes only five inputs (not MobiDB, InterPro, CD-CODE, the history or the CPLM directory). No R `sessionInfo()` is saved.
7. **The deployed website data has no provenance.** `docs/data/methylarg_sites.json` is the published union, exported 2026-10-16 and patched in place since (CPLM, CARM1, HPA, disease families, Maron calls) without updating `generated`. It has no `meta.provenance` and no reference checksums, so the explorer's sequence-drift guard never fires. `methylarg_evidence.json` has no metadata.
8. **The site mixes reference versions and dates.** Site positions follow the April snapshot; `scripts/37` checked crosslinks, IntAct and control arginines against the UniProt it downloaded at export time; the browser draws live UniProt, AlphaFold and PDB numbering. Condensate membership on the site (CD-CODE export of 2026-03-18) and the per-protein CD-CODE details (API release 2.3; now 2.4) disagree for 5 condensates. Disorder classes come from MobiDB; the maps draw live UniProt disordered regions.
9. **Maron Table S2/S3 are not staged anywhere,** so `scripts/16` cannot be rerun from the repository.

## 5. Clean-up proposal (for approval; nothing done yet)

Goal: a repository anyone can rerun, with results unchanged. Each step is checked by two tests that must keep passing: (a) the pre-fix remapping on `inputs/` reproduces the published remapping on 28,960 of 28,960 rows; (b) the corrected pipeline on `inputs/` gives the same union and the same downstream tables as the run being kept.

1. **Inputs.** Point `hpc/ptm_env.sh` at `inputs/` (the files now have their final names), add a step that refuses to run when a file's checksum differs from `inputs/MANIFEST.tsv`, and let the hygiene guard accept exactly the manifest's files. Stage working copies elsewhere so no step writes into `inputs/`.
2. **Complete the frozen set.** Rebuild the InterPro intervals for all reviewed proteins from the April `protein2ipr.dat.gz` and add them; the accession history and Maron Tables S2/S3 are now in `inputs/`. Put the large raw files (InterPro 16 GB, MobiDB JSON, CPLM, the `scripts/37` cache) in a Zenodo record and list them in the manifest by checksum.
3. **Fix the hazards** above: no empty InterPro files from step 00; explicit inputs for script 24; a run record (input checksums, commit, R `sessionInfo()`) written by every run; pin the R data packages.
4. **Layout.** One provenance document (this file, `DATA_PROVENANCE_REVIEW.md` and `ISSUES_2026-10-19.md` merged), the README reduced to how to run; retire the legacy `data/` and `results/` stubs; keep `PTM_results/` as the published release, labelled as such.
5. **Then** the rerun and the website rebuild (`ISSUES_2026-10-19.md` issue 17), the other-PTM export fix, and a pull request to `main`.

Decisions for the lab: issue 9 (paralog windows, 22 rows), whether the website should show isoform-numbered other-PTM sites after a window check, whether to extend the master beyond iPTMnet (changes cross-PTM results), and which release the website serves by default.

## 6. Effect on the paper's featured proteins (release v2.0 vs the published union, v1.0)

| Protein | Published sites | Release sites | Removed | Added |
|---|---|---|---|---|
| FUS (P35637) | 30 | 34 | — | R422, R518, R522, R524 (peptide recentred) |
| FMR1 (Q06787) | 13 | 13 | R562, R564 (evidence mismatch: windows from an isoform-6-only exon) | R483, R558 (peptide recentred) |
| Histone H3.1 (P68431) | 10 | 10 | — | — (R64 flagged `paralog_window`: its window matches H3.1t exactly) |
| SmD1 / SNRPD1 (P62314) | 6 | 14 | — | R92, R104, R108, R110, R112, R114, R118, R119 (peptide recentred: the C-terminal RG tail) |
| SmD3 / SNRPD3 (P62318) | 6 | 7 | — | R125 (peptide recentred) |
| SmB / SNRPB (P14678) | 13 | 13 | — | — |
| EWSR1 (Q01844) | 35 | 36 | — | R518 (peptide recentred) |

None of the 117 sites removed by peptide recentring belongs to these proteins. The published architecture figures under-draw SmD1 most: the published union held 6 of its 14 sites.

## 7. The other session's issues and README, checked on the April snapshot (data and provenance only)

`ISSUES_2026-10-19.md` (counted there on UniProt 2026_03):

| Issue | On the April snapshot |
|---|---|
| 18. Records on retired accessions whose evidence is in no current entry | Confirmed: 315 unplaced records (Larsen 176, Musiani 106, PhosphoSitePlus 13, Lim 10, Maron 4, Fedoriw 2), all `unresolved_no_anchor`. The UniSave fix would place some of them and change the union. |
| 19. Deleted accession, one paralog holds the evidence | Confirmed: TMEM189 G3V2F7, CALM2 H0Y7A7, ELMO1 Q5JVZ5 and SMARCA4 Q9HBD4 each have unplaced records. A policy question with issue 9. |
| 20. Evidence check laxer than step 02 for proteome-wide placements | A future risk, not a present error: all 15 proteome-wide placements agree with the check. |
| 21. Accession history from a later release than the FASTA | Confirmed and recorded in `inputs/MANIFEST.tsv` (2026_03). Its effect is 3 records: with it HIST2H2BE U3KQK0 and CTAGE5 O15320 are placed and the TMEM189 readthrough row (G3V2F7) is not; without it, the reverse. |
| 22. Identical paralogs from split accessions | Confirmed (HSPA1A/HSPA1B, the TUBA3C pair, the PRAMEF15 pair are identical in the snapshot). More generally, 593 release sites (4.4%) sit in a 31-residue window that occurs identically in another reviewed entry; the release flags them `shared_sequence` and names the other entries. |

README statements about data and provenance:

- *"The job reuses the staged inputs … so no newer UniProt or InterPro release enters the comparison."* Not so. The first step empties the staged InterPro intervals (reproduced with `main`'s code and its pinned pandas), so step 03 stops unless the intervals are read from elsewhere; the accession history is downloaded from UniProt `current_release` when missing; and the intervals are restaged from `protein2ipr.dat.gz`, or from InterPro's current release if that file is absent.
- *"Use `hpc/slurm_run_integration.sh` as the authoritative ordered recipe."* The published results did not come from it with its defaults: both clustering steps used 100 permutations, not 100,000 (`clustering_qc.tsv`, the null-permutation tables); step 03 read InterPro intervals that step 00 had not emptied; and the tables were written by pandas 3, not the pinned 2.2 (blank cells written as empty rather than `nan`). `main`'s environment file cannot be installed as written (its Bioconductor packages need the bioconda channel).
- *"DOSE 4.4.0 with org.Hs.eg.db 3.22.0 reproduces the committed family tables exactly."* Being checked: the cloud rerun uses these versions.
- *CPLM 4.0 validated against UniProt release 2026_03.* The lysine sites are checked against a different reference than the methylarginine union (the April snapshot), so the lysine-methylation layer and the K/R co-occurrence analysis (`scripts/40`) mix two references; the 13 proteins whose sequence differs between them are the exposure. Validate CPLM against the same FASTA.
