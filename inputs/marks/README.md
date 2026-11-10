# Per-mark inputs

One staged table per mark and source, for the context marks of the release (`releases/v1.5.0/`). Columns: `canonical_UniProtAC`,
`position`, `residue`, `mark`, `state`, `evidence`, `evidence_class`, `pmids`, `source`, `source_release`. Protein
keys and positions are on the reference proteome (`inputs/reference/`, UniProt 2026_03).

| File | Holds | Staged with |
|---|---|---|
| `uniprot_<release>_methyllysine_sites.tsv` | N6-methyllysine (me1), N6,N6-dimethyllysine (me2), N6,N6,N6-trimethyllysine (me3), N6-methylated lysine (me: degree not stated) | `scripts/54_stage_uniprot_ptm_sites.py` |
| `uniprot_<release>_phosphorylation_sites.tsv` | phosphoserine, phosphothreonine, phosphotyrosine (state `phospho`) | the same |
| `uniprot_<release>_methylarginine_crosscheck.tsv` | UniProt's arginine methylation, kept only to cross-check the methylarginine release | the same |
| `iptmnet_<release>_methyllysine_sites.tsv`, `..._phosphorylation_sites.tsv`, `..._methylarginine_crosscheck.tsv` | iPTMnet 6.2 human records (CC BY-NC-SA 4.0), degree not stated (state `me`); `upstream_sources` keeps the upstream database of each record | `scripts/55_stage_iptmnet_ptm_sites.py` |
| `cplm_<release>_methyllysine_sites.tsv` | CPLM 4.0 human lysine methylation sites that match the reference sequence (state `me`: no degree), with CPLM's evidence label (`cplm_evidence`: Exp., Exp./Dat., Dat.) and its references verbatim; `evidence_class` experimental for Exp. and Exp./Dat., unspecified for Dat.; `pmids` the references that are PubMed IDs | `scripts/39_stage_cplm4_lys_methylation.py --marks-dir inputs/marks` |
| `berryhill2023_methyllysine_validated.tsv`, `..._rejected.tsv`, `..._summary.json` | Berryhill et al. 2023 human lysine methylation sites (me1, me2, me3; HEK293T, U2OS), checked on the reference proteome; see the last section for the licence note | `scripts/57_stage_berryhill_lys_methylation.py` |
| `uniprot_ptm_sites_summary.json`, `iptmnet_ptm_sites_summary.json`, `cplm_ptm_sites_summary.json` | release served, rows, drops, MD5 of the input and of each table | the same |
| `raw_files.tsv` | URL, size, MD5, SHA-256, date of each raw download (kept in `PTM_data/marks/`, not in Git) | |

UniProt is CC BY 4.0. iPTMnet's license page says CC BY-NC-SA 4.0 while its download readme also mentions CC BY 4.0; the stricter terms are followed, as the release is CC BY-NC-SA 4.0. Its API serves the current release only (2026_04 on 2026-11-09), newer than the reference
(2026_03): an entry whose sequence differs from the reference is left out whole, so no position moves between
sequences (none differed on this run).

## CPLM 4.0

Zhang et al., Nucleic Acids Res 50:D451-D459 (2022), doi:10.1093/nar/gkab849; `Methylation.zip` from cplm.biocuckoo.cn (its `Methylation.txt` is dated 2021-06-13; URL, date and checksums in `raw_files.tsv`). **CPLM states no license**: the site carries a copyright notice, "All Rights Reserved". On 2026-11-09 (confirmed 2026-11-10) the owner decided to include its human lysine methylation sites in the release, cited, pending the authors' permission; the terms the authors give are followed. The download itself is not redistributed. Rows are kept when the lysine is in CPLM's own sequence and in the reference sequence at the same position, and the sequences are identical or the +/-7 residue windows agree; nothing is remapped (rejected rows with their reasons are counted in `cplm_ptm_sites_summary.json`). Lysine acylations of CPLM are not staged.

## Berryhill et al. 2023: staged

Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. Global lysine methylome
profiling using systematically characterized affinity reagents. Sci Rep 13:377 (2023),
doi:10.1038/s41598-022-27175-x (PMC9825382). **Licence: article CC BY 4.0; supplementary files carry no separate notice;
included by the owner's decision (2026-11-10), without waiting for the authors' confirmation.** Cite the paper wherever
these sites are described.

`scripts/57_stage_berryhill_lys_methylation.py` reads Sup Table 5 of Supplementary Information 2
(`..._MOESM2_ESM.xlsx`, listed with URL, date and SHA-256 in `raw_files.tsv`; the download stays out of Git) and writes
`berryhill2023_methyllysine_validated.tsv`, `berryhill2023_methyllysine_rejected.tsv` and
`berryhill2023_methyllysine_summary.json` (counts and checksums). Not used: the other sheets (Sup Tables 4 and 7 carry
PhosphoSitePlus information) and the other columns of Sup Table 5 (PSM counts, motif, description) and Supplementary
Information 3.

A row is kept when the canonical sequence has K at the position and the peptide lies there; otherwise the peptide is
looked up in the canonical sequence and the site moves to its K when the match is unique (`relocated`, with
`source_position`); otherwise the row is rejected with a reason. States m1, m2, m3 are me1, me2, me3. Columns are the ten
of the other mark tables (`evidence_class` is `high-throughput`, `source` is `Berryhill2023`) plus
`site_probability_max`, `cell_lines`, `enrichments`, `peptides`, `relocated`, `source_position`. The paper filters peptides
for "high confidence identification and modification site probability" but gives no number in the text; the column is
kept, not filtered (the table's values are 90.1-100; PTMScan Rep 2 rows have none, so their sites may have an empty value).

The validated table is not named `*_sites.tsv` on purpose: `scripts/56_build_context_marks.py` globs that pattern and has
no `SOURCES` entry for Berryhill2023, so the merge into a release is a later, deliberate step.

## Ochoa et al. 2020: staged, not yet in a release

Nat Biotechnol 38:365-373, doi:10.1038/s41587-019-0344-3. Added by the owner's decision (2026-11-10). The article is
not open access: its page says "© 2019 The Author(s), under exclusive licence to Springer Nature America, Inc.", and
the supplements carry no separate notice. The staged tables are therefore `checksum_only` (not in Git, listed in
`.gitignore`); `scripts/58_stage_ochoa_phosphosites.py` rebuilds them from the downloads (Supplementary Table 3,
`..._MOESM5_ESM.xlsx`, for the functional score; Supplementary Table 2, `..._MOESM4_ESM.xlsx`, for the residue and
the evidence columns). The counts file is in Git.

| File | Holds |
|---|---|
| `ochoa2020_phosphorylation_validated.tsv` | phosphoserine, phosphothreonine and phosphotyrosine sites on the reference proteome, one row per residue: the ten columns above (source `Ochoa2020`, evidence class `high-throughput`) plus `functional_score`, `best_localization_prob`, `best_PEP`, `biological_samples`, `spectral_counts` |
| `ochoa2020_phosphorylation_rejected.tsv` | sites that failed the check, with `reason` and `detail` |
| `ochoa2020_phosphorylation_summary.json` | counts, drops by reason, score quantiles, MD5 and SHA-256 of the inputs |

A site is kept when its accession is a canonical entry of the reference and the sequence has the table's residue at
the position. Neither table gives a peptide or window, so nothing is relocated: a mismatch is rejected. The file
names end in `_validated`, not `_sites`, so `scripts/56` does not read them until the release merge adds them.
