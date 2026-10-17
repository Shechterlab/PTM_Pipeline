# Data provenance review: methylarginine site union

Handover notes for a fresh review of how source records become union sites. Written 2026-10-17 on branch `claude/awesome-thompson-7j9n92`. Everything below was checked against the raw source files and current UniProt; the numbers can be reproduced with the scripts named in each section.

## Where things are

| What | Where |
|---|---|
| Raw sources | Not in the working tree; in git history at commit `cd56dd0`: `data/base/human_ptm_master.tsv`, `data/external/Table_S5_Compiled_Methylarginine_Data.xlsx` (Maron et al. 2021 Table S5), `data/external/mmc2.xlsx` (ProMetheusDB). Extract with `git show cd56dd0:<path> > file`. On the HPC they are staged under `PTM_DATA_ROOT`. |
| Parsers | `scripts/01_parse_maron_s5.py`, `scripts/01_parse_prometheus_mmc2.py` (master rows are read directly in step 02) |
| Remap to canonical UniProt | `scripts/02_remap_arg_methyl_to_canonical.py` |
| Union and deduplication | `scripts/02_integrate_arg_methyl_sources.py` |
| Published results | `PTM_results/` (13,719 sites on 5,090 accessions), used in DeAngelo, Silverstein & Shechter, *Cell Chem Biol* 33:1404–1421 (2026), [doi:10.1016/j.chembiol.2026.09.011](https://doi.org/10.1016/j.chembiol.2026.09.011) |
| Website data | `docs/data/*.json`, exported from `PTM_results/` by `scripts/36_export_site_explorer_data.py` (and `37_` for interactions). **Still the published union.** |
| Rerun after the fixes | `sbatch hpc/slurm_rerun_remap_fix.sh /path/to/PTM_Pipeline` → new results root, `rerun_comparison/report.md`, evidence checks of both roots |
| Checks | `scripts/qc_remap_evidence_check.py` (every row against its own evidence, independent of step 02), `scripts/qc_remap_hand_check.py` (one protein by eye; `--audit-union` residue audit), `scripts/38_compare_results_roots.py` (table-by-table diff of two results roots) |

## What is already verified

- **Parsers are faithful.** Parsed rows match the raw files record for record: Maron S5 26,515 of 26,519 (the 4 dropped rows have no accession), ProMetheusDB 1,598 of 1,598 R rows, master table 847 of 847.
- **Maron S5 reproduces the Maron paper.** 15,381 unique (accession, position) on 5,252 accessions (paper: 15,386 on 5,255; the gap is the 4 rows without accession). Maron's own data: 1,829 sites on 585 proteins, of which 630 sites and 106 proteins are not in any other S5 study, exactly as the paper states ([Maron et al. 2021, iScience 24:102971](https://doi.org/10.1016/j.isci.2021.102971)).
- **The published run's sequences equal today's UniProt** at 12,845 of 12,847 sites with a saved window (`PTM_results/motif_arg_odds/sequence_windows.tsv`). Exceptions: RAB12 (Q6IQ22; UniProt extended its N-terminus on 2026-09-02, after the run) and a "site" at STAU1 (O95793) R946, beyond the protein's end (see issue 6).
- **Replaying the old step 02 on today's sequences reproduces the published remapping** on 28,958 of 28,960 rows (the two others are RAB12).
- **Sequence drift since the sources is real but small.** 26,204 of 28,960 rows (90%) have a reported position that still fits the current sequence exactly; 1,565 have their evidence elsewhere in the same entry, 674 have an accession that no longer exists, 505 have evidence not in that entry, 12 sit on a non-R residue. Most drift is Larsen 2016. 297 source accessions no longer exist (159 merged, 138 deleted).

## Issues

Severity: **high** changes sites in the published union; **medium** loses or keeps doubtful sites; **low** bookkeeping. "Fixed" means fixed in code on this branch; nothing in `PTM_results/` or `docs/data/` has been regenerated yet.

| # | Issue | Rows / sites | Status |
|---|---|---|---|
| 1 | Peptides read as centred windows: the reported arginine was replaced by the one nearest the peptide middle | 1,773 rows | high · fixed |
| 2 | Position computed on the reviewed entry but filed under the source's unreviewed accession | 857 rows; 578 sites counted twice | high · fixed |
| 3 | Lysine-centred Musiani windows turned into arginine sites | 131 rows (268 K windows in S5) | high · fixed |
| 4 | Windows cut off at a protein end centred wrongly | 302 rows | high · fixed |
| 5 | Guessed arginine for multi-R peptides whose reported position does not fit | 52 rows | medium · fixed (no longer placed) |
| 6 | Filed under the wrong reviewed protein | 49 rows | high · fixed |
| 7 | Candidate search by gene symbol misses outdated symbols and deleted accessions | 129 rows dropped | medium · **open** |
| 8 | Reported arginine kept although the source window does not match there | 139 rows | medium · **open** |
| 9 | Larsen windows that match a paralog, not the reported protein | 31 rows | medium · **open** |
| 10 | Malformed source accessions | 9 rows | low · **open** |
| 11 | No output validation in the pipeline | — | process · partly addressed |
| 12 | Downstream statistics bugs found along the way | — | fixed |

Totals with the fixes (local rebuild, same sequences): 13,354 sites on 4,712 reviewed proteins; 12,512 sites kept, 842 added, 1,207 removed; 80 reviewed proteins enter, 28 leave. Every corrected site is an arginine on a reviewed sequence (published: 872 sites without a usable sequence and 21 on non-arginines). The independent check agrees with the corrected placement on 27,677 rows and contradicts it on 31 (issue 9); for the published run it agrees on 24,757 and contradicts 3,028.

### 1. Peptides read as centred windows (fixed: commits `1957827`, `2b13e0a`)

Maron S5's "Peptide Sequence" column holds centred windows for six studies (Guo 13 aa; PhosphoSitePlus, Musiani, Fong 15; Larsen, Wei 31) and whole peptides for four (Maron's own data, Fedoriw, Li, Lim). Step 02 treated every token as a window: the middle residue was the site, or, if it was not R, the R nearest the middle (`common.centered_anchor_indices`). For peptides this ignored the reported position, even when that position was an R inside the peptide in the same sequence. In all 1,773 affected rows the reported position fits today's sequence, and in 1,767 the protein sequence has not changed since before the source was published, so this is not drift. Multi-site peptides collapsed onto one site (SNRPD1 R104–R119 onto R98–R106; HNRNPH1 R206, R212, R224 onto R217; TAF15 R467, R483, R490 onto R475). The published union holds only 288 of Maron's 630 new sites at their reported positions; the corrected code keeps 620.

Fix now in `scripts/02`: the kind of token is read per source from the data (`annotate_window_lengths`: a source whose tokens mostly share one window length gives windows); peptides keep their reported arginine (`sequence_position_verified`) and place a site alone only when they contain a single R. A first version decided by token length (13/15/31) and broke 39 peptide rows; that is fixed in `2b13e0a`.

*For review:* confirm the window/peptide classification per reference (printed by step 02) and that no other source mixes the two.

### 2. Filed under the unreviewed accession (fixed: `0a9cf10`)

When a row's own accession was unreviewed (TrEMBL), the gene-name candidate search found the peptide on the reviewed entry and computed the position there, but `corrected_canonical_UniProtAC` stayed the row's accession. Result: reviewed-entry numbering on unreviewed accessions (A2A274 R474 is a G in A2A274; ACO2 Q99798 R474 is the site), no sequence for motif/disorder analyses (`motif_analysis_qc.tsv` reports 872 sites "without usable sequence"), and 578 sites duplicated across the reviewed and unreviewed accession. Fix: file under the matched accession (`is_accession_remapped`).

### 3. Lysine-centred windows (fixed)

268 of Musiani's 1,712 rows in Maron S5 are 15-mer windows centred on a lysine, with a lysine at the reported position (sequences unchanged since 1998–2007). Musiani et al. 2019 measured both R and K methylation ([doi:10.1126/scisignal.aat8388](https://doi.org/10.1126/scisignal.aat8388)); these look like lysine sites carried into an arginine compilation. The old code turned 131 of them into an R a few residues away. ProMetheusDB (same lab) lists 11 of the 131 as methyl-lysine sites.

*For review:* check Musiani's original supplementary table to confirm these are lysine sites, and consider logging them at the parse step rather than silently giving no anchor.

### 4. Truncated windows (fixed: `2b13e0a`)

Windows shorter than full length (cut at a protein end, or at the ubiquitin C-terminus inside RPS27A/UBA52) were centred on their own middle. The site is half a window from the intact end; both readings are now tried and only an arginine counts (`window_anchor_indices`).

### 5. Multi-R peptides that do not fit (fixed)

52 rows (Fedoriw, Lim, Maron, Li) whose peptide is in the protein but does not cover the reported position, and contains several arginines. Nothing in the evidence says which arginine; the old code took the middle one. Now unplaced.

### 6. Wrong protein (fixed: `0a9cf10`)

Same mechanism as issue 2 with a reviewed accession, for example a PTCD1 window filed as MED12 (Q93074) "R235" (a Y in MED12), and a Lim XRN2 peptide filed under STAU1 (O95793) at 946, beyond STAU1's 577 residues (now XRN2 Q9H0D6 R946).

### 7. Recall gap from gene-symbol candidates (open)

`candidate_accessions` searches the row's own accession plus accessions sharing its gene symbol. Rows with an outdated symbol (MST4 → STK26, FAM60A → SINHCAF, C11orf67 → AAMDC) and an unreviewed or deleted accession get no candidates and are dropped, in both versions. A proteome-wide search finds 129 of them on a single reviewed entry.

*Suggested fix:* when no candidate matches, search full-length windows and marked peptides across the whole reviewed proteome and accept a unique hit (`qc_remap_evidence_check.Proteome` shows the approach). Map outdated symbols through UniProt synonyms as well.

### 8. Fallback keeps a reported arginine the evidence does not support (open)

`fallback_position` keeps the reported position whenever that residue is an R in the reported entry, even when the source window is not found there. Of 139 such rows from window sources, the window differs from the sequence at the reported position by 0–2 residues in 54 (variants or conflicts; probably right), 3–6 in 46, and more than 6 in 39 (35 of them on sequences unchanged since the source, so the source paired the window with the wrong protein or position).

*Suggested fix:* for window sources, require the window at the reported position to match with at most 2 mismatches; otherwise leave the row unplaced or flag it.

### 9. Paralog windows (open)

31 Larsen rows (and 1 Fedoriw) whose 31-residue window matches a paralog exactly (APP vs APLP2, PFKM vs PFKL, HSP90AB1 vs HSP90AB3P) and not the reported protein, while the reported position is an R. Both versions keep the reported accession. This is the shared-peptide problem; the window is the stronger evidence for which sequence was observed.

*Suggested fix:* when the window matches another reviewed entry exactly and not the reported one, place the site there, or keep it and flag it as ambiguous.

### 10. Malformed accessions (open)

`P31943; P55795` (4 rows, two accessions in one cell; placed on P31943), `NP_001076783` (4 rows, RefSeq; unresolved), `Q9UFP1_VAR_Q460R` (1 row, variant suffix; unresolved). The 4 S5 rows with no accession are also dropped; their windows could place them.

### 11. No output validation (partly addressed)

QC counted coverage, not correctness. Signals that were printed but read as expected drift: "rows_with_position_shift 3,954" (14% of rows) in `remapped/summary.tsv`, and 872 sites "without usable sequence" in the motif QC.

*Suggested fix:* make `scripts/qc_remap_evidence_check.py` and `qc_remap_hand_check.py --audit-union` standard pipeline steps that fail the run when any union site is not an R, sits on an accession outside the reviewed FASTA, or contradicts its evidence beyond an agreed list. The rerun script already runs the evidence check.

### 12. Downstream fixes found along the way (fixed)

- `scripts/08` grouped per-protein arginine odds by accession and gene symbol, splitting 138 proteins with alias symbols (MRE11/MRE11A) into weaker tests; q ≤ 0.05 proteins go from 807 to 832 on the published union (`0a9cf10`).
- `scripts/08` drew heatmap stars from raw p while the Fig 3A legend states BH q; on the published data they differ at one cell, P +1 (`5b4a272`).
- The repo hygiene guard fails on a checkout that tracks `PTM_results/`; `PTM_SKIP_HYGIENE=1` skips it.

## Effect on the paper (published vs corrected, local rebuild)

No conclusion changes. Fig 3A: 12 of 84 cells change significance class; no significant cell changes direction. Fig 3B disordered log2 OR 1.21 → 1.46, ordered −1.39 → −1.29, boundary −0.26 → −0.20, "no annotation" 877 → 5 sites. Fig 3C (GO) and 3D unchanged in substance (FUS gains R422, R518, R522, R524; FMR1 gains R483, R558). Motif families: RGG 4.17 → 4.35, RG 2.55 → 2.63, DR 1.49 → 1.48, GAR 0.59 → 0.71 (log2 OR vs arginines in methylated proteins). Independent of the bug: "most of these sites are found within IDRs" holds only counting IDR edges (strictly inside: 45% published, 50% corrected), and the 3B legend mentions protein numbers the panel does not show.

## Suggested order of work

1. Review issues 7–10 and decide the policy for each; implement in `scripts/02` with the evidence check as the test.
2. Run `hpc/slurm_rerun_remap_fix.sh` on the staged inputs; read `rerun_comparison/report.md` and the two evidence-check reports.
3. Regenerate the website data (`scripts/36_export_site_explorer_data.py`, then `37_` with `PTM_RERUN_INTERACTIONS=1`) and the Zenodo archive (new version with a changelog citing the comparison report).
4. Make the evidence check and residue audit permanent pipeline steps (issue 11).
