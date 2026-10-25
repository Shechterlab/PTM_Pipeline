# Changes from v1.3.0 to v1.3.1

| | v1.3.0 | v1.3.1 |
|---|---:|---:|
| Sites | 11,601 | 11,601 |
| Proteins | 4,432 | 4,432 |
| Sites kept pending | 0 | 0 |
| Source records | 66,660 | 66,660 |
| Source records placed | 37,518 | 37,518 |
| Studies and databases | 18 | 18 |

Of the 11,601 sites of v1.3.0 (pending ones included), 11,601 are kept and 0 removed; 0 sites were added. The records of v1.3.0 are read in the same order and all are placed as before. 

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## Corrections after the review of v1.3.0

v1.3.1 keeps the sites, records and placements of v1.3.0. It changes how some evidence is called and described, after an outside review of v1.3.0.

### Inhibitor and knockout responses: one rule, read at the peptide level

Maron et al. 2021 Table S3 and the CARM1 knockout of Shishkova et al. 2017 measure methylated peptides, not sites. A site's call over the peptides that carry it is now `decreased` when at least one peptide decreased and none increased, `increased` likewise, `mixed` when peptides changed in opposite directions, and `no_significant_change` (in place of `no_change`) when none changed: measured, not shown to be unchanged. The release, the website and the Maron 2021 summaries of the analyses (`scripts/16`) now share this rule (`scripts/peptide_calls.py`); the summaries had called a site from its median change, which can have the opposite sign to a significant peptide. The median stays as a description only. Every site keeps its v1.3.0 call; only `no_change` is renamed.

| Sites called | GSK591 (PRMT5 inhibitor): v1.3.0 → this version | MS023 (type I PRMT inhibitor): v1.3.0 → this version | CARM1 knockout: v1.3.0 → this version |
|---|---:|---:|---:|
| `decreased` | 159 → 159 | 196 → 196 | 257 → 257 |
| `increased` | 92 → 92 | 177 → 177 | 2 → 2 |
| `mixed` | 46 → 46 | 58 → 58 | 1 → 1 |
| `no_significant_change` | 626 → 626 | 492 → 492 | 118 → 118 |

- A peptide that carries several methylarginines reports them together, so its change cannot be given to one of them. The tables now count those peptides (`n_multi_site_peptides`) and give the call from the single-site peptides alone (`call_single_site_peptides`, empty when there are none). Of the 1,049 site, antibody and inhibitor rows of `supporting/maron2021_inhibitors` with a change, 884 rest only on peptides carrying another methylarginine (no single-site peptide, or none of those changed).
- Most rows of Table S3 hold imputed values; where one condition is a constant imputed value, the peptide was detected in the other condition only. `n_imputed_peptides` and `n_one_condition_peptides` count them; 175 of the rows with a change have every peptide detected in one condition only.
- CARM1 knockout: a peptide now changes at Benjamini-Hochberg q < 0.05 (the study's p values adjusted over the quantified peptides of each cell line) and a twofold change; the study's own rule (nominal p < 0.05) is kept as `call_source_rule`. Both rules call the same 363 site and cell-line rows changed, each the same way.
- Neither study measured total protein in the same samples, so a peptide's change can also be a change in its protein. The website calls a decrease after a PRMT inhibitor "Writer suggested", not a writer.

### RBR-ID: called on the biological samples

RBR-ID (Oksuz et al. 2023) measured five biological samples per condition, each in up to three runs. v1.3.0 followed the published rule, which counts every run as an independent observation and calls a peptide at unadjusted p < 0.05. The call is now made on the samples: each peptide's values are fitted for sample and acquisition date, compared with Student's t test (at least three samples per condition), and called at Benjamini-Hochberg q < 0.05 with a decrease on crosslinking. The published rule is kept beside it as a reproduction (columns ending in `_published` in `inputs/interactions/`, and the test category `rbrid_residue_published_rule`).

| | Biological samples, BH (used) | Published rule |
|---|---:|---:|
| RBR-ID+ peptides | 3,039 | 9,230 |
| Residues under them | 41,971 | 121,523 |
| Methylarginines against the other arginines of the same proteins (Mantel-Haenszel) | OR 2.09, p 1.5e-21, 488 proteins | OR 1.95, p 8.4e-40, 925 proteins |

The residue score threshold of the RBR-ID regions is now set by a rule fixed in advance: the lowest of 1, 2, 3 and 5 whose permutation false-discovery rate is under 10 % (2; v1.3.0 used 3).

With fewer RBR-ID+ peptides, 462 sites lose `rbrid_residue` from `interaction_evidence` and 5 gain it. 214 sites change `interaction_status`: 184 from experimental to none, 27 from experimental to predicted_only, 3 from none to experimental. These are the sites whose only experimental evidence of RNA or DNA contact was, or now is, RBR-ID.

### How each site is localised

Sites rest on different kinds of evidence, which do not share one error rate. The site table now names them: `localised_by` (a localisation score, a peptide that leaves one possible arginine, or a curated annotation) and `curated_only`, true for the 1,314 of 11,601 sites that only a curated annotation (UniProt with experimental evidence, PhosphoSitePlus as compiled in 2021) places. The website can leave them out (Filters, Provenance).

### DisGeNET terms tested with clusterProfiler

The DisGeNET terms are now tested like every other protein or gene set, with clusterProfiler's enricher, on DOSE 4.4.0's DisGeNET term-gene table (`inputs/disease/dose_4.4.0_disgenet_term_genes.tsv`). That is the test DOSE's enrichDGN ran on the same data, which later DOSE releases removed. Every count, p value and q value of the 8,334 tested terms is unchanged.

### Release files and the website

- The checksums are written after every other file, and the build checks them (`MD5SUMS` here and in each analyses folder); in v1.3.0 the two analyses READMEs were rewritten after their checksums.
- `release.json` records the Python and package versions that wrote the release (`python_environment`).
- Website: a test run on a selection made by the feature it tests is shown as descriptive and left out of the multiple-testing correction; other selection tests are labelled exploratory. Proteins without a site get their 3D structures and, when they are in the interaction data, their STRING partner network.

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (Site specificity in the CHANGES.md of v1.2.0): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.3.0_to_v1.3.1.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.3.1.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.
