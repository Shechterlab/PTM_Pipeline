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
| `uniprot_ptm_sites_summary.json`, `iptmnet_ptm_sites_summary.json` | release served, rows, drops, MD5 of the input and of each table | the same |
| `raw_files.tsv` | URL, size, MD5, SHA-256, date of each raw download (kept in `PTM_data/marks/`, not in Git) | |

UniProt is CC BY 4.0. iPTMnet's license page says CC BY-NC-SA 4.0 while its download readme also mentions CC BY 4.0; the stricter terms are followed, as the release is CC BY-NC-SA 4.0. Its API serves the current release only (2026_04 on 2026-11-09), newer than the reference
(2026_03): an entry whose sequence differs from the reference is left out whole, so no position moves between
sequences (none differed on this run).

## Berryhill et al. 2023: downloaded, not staged

Sci Rep 13:377, doi:10.1038/s41598-022-27175-x (PMC9825382). The article is CC BY 4.0. The two supplements
(`..._MOESM2_ESM.xlsx`, `..._MOESM3_ESM.xlsx`) have no license notice of their own, the journal's license policy
speaks of "the article", and one sheet (Sup Table 7) mixes the authors' sites with PhosphoSitePlus sites, so the
supplements' terms are not confirmed. They are listed in `raw_files.tsv` and kept out of Git; no table is staged
from them. Sup Table 5 (6,576 rows, 5,089 sites on 2,751 proteins) gives protein, position and state (m1, m2, m3),
and would be staged as a second methyllysine source (`source=Berryhill2023`, `evidence_class=high-throughput`) once the
terms are confirmed in writing by the authors or Springer Nature. Do not use Sup Tables 4 and 7: they carry
PhosphoSitePlus information.
