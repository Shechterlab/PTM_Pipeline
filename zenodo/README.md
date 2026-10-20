# PTM_Pipeline v1.1.0: offline and Zenodo archive

This archive holds the data behind the Methylarginine Site Explorer
(<https://shechterlab.github.io/PTM_Pipeline/>) and the inputs of PTM_Pipeline v1.1.0. It was made on 20 October 2026
from branch `claude/awesome-thompson-7j9n92` of <https://github.com/Shechterlab/PTM_Pipeline>.

It covers:
- the readable tables;
- every raw download they were made from;
- a complete copy of the Git repository.

Please cite DeAngelo JD, Silverstein AM, Shechter D. *Protein arginine methylation rewires molecular interactions.*
Cell Chem Biol 33(10):1404–1421 (2026), doi:10.1016/j.chembiol.2026.09.011, and the original source of any data you
use.

## Files

| File | Size | Destination | Contents |
|---|---:|---|---|
| `ptm_pipeline_inputs_tables_2026-10-20.tar.gz` | 45.0 MB | Public | inputs/reference, inputs/interactions, inputs/layers as in Git |
| `ptm_pipeline_interactions_raw_cache_2026-10-16.tar.gz` | 160.7 MB | Public | The raw download cache of scripts/37: 36,550 files |
| `ptm_pipeline_layers_raw_downloads_2026-10-18.tar.gz` | 15.5 MB | Public | HPA 25.1, Shishkova 2017 Supplementary Data 1, the UniProt 2026_03 FASTA |
| `ptm_pipeline_RESTRICTED_cplm4_2026-10-18.tar.gz` | 2.5 MB | Restricted | The CPLM 4.0 download and its tables (all rights reserved) |
| `PTM_Pipeline_all_branches_2026-10-20.bundle` | 176.5 MB | Offline | The Git repository: 6 branches, 3 tags, full history |
| `README.md` |  | Public | This file |
| `zenodo_public.json, zenodo_restricted.json` |  | Upload | Zenodo deposit metadata |
| `zenodo_upload.py` |  | Upload | Creates a Zenodo draft, uploads and checks the files |
| `github_release.zenodo.json` |  | Upload | Metadata for the GitHub v1.1.0 release (save as .zenodo.json) |

**Where each file goes:**
- **Public:** a public Zenodo record (`zenodo_public.json`).
- **Restricted:** a Zenodo record with restricted access (`zenodo_restricted.json`), or private storage only.
- **Offline:** lab storage. Everything in the bundle is already public on GitHub, but the bundle is the one file
  that restores the whole repository, history included.

`MD5SUMS` and `SHA256SUMS` cover every file above.

## Check the files

```
sha256sum -c SHA256SUMS          # Linux
shasum -a 256 -c SHA256SUMS      # macOS
```

## What each archive holds

- **`ptm_pipeline_inputs_tables_2026-10-20.tar.gz`** unpacks to `inputs/reference/`, `inputs/interactions/` and
  `inputs/layers/`, exactly as in Git at commit `b75633e`. Each folder has a README and a `MANIFEST.tsv` (rows,
  MD5, source URL, version, licence, how the table was made).
  - **reference:** UniProt release 2026_03, the reviewed human reference proteome UP000005640 (20,416 sequences).
    It is the one UniProt reference for the pipeline, with its entries table and accession history.
  - **interactions:** 28 readable tables and the UniProt FASTA the export used, from these sources:
    - UniProt, QuickGO, STRING v12.5, Complex Portal;
    - PDBe-KB, SIFTS, PDB contacts and interface geometry;
    - ProtVar, AlphaFold DB, IntAct;
    - pepR-MS, RBDmap, Trendel 2025, CD-CODE 2.3.

    It also holds `raw_cache_files.tsv` (the MD5 of each of the 36,550 cached files).
  - **layers:** HPA 25.1 subcellular locations, the CARM1 knockout response (Shishkova et al. 2017), and the CPLM
    summary. The CPLM tables themselves are in the restricted archive.
- **`ptm_pipeline_interactions_raw_cache_2026-10-16.tar.gz`** unpacks to `interactions_cache/`, the download cache
  of `scripts/37_export_interaction_data.py`, saved 16 to 19 October 2026. It holds:
  - bulk files from UniProt, QuickGO, STRING, Complex Portal, SIFTS, ProtVar and IntAct;
  - one API answer per protein or structure, from PDBe-KB, AlphaFold DB, RCSB, ProtVar and CD-CODE;
  - the 10 Å atom neighbourhoods measured for interface geometry;
  - the pepR-MS, RBDmap and Trendel 2025 spreadsheets.

  With this cache in place and network access off, `scripts/37` rebuilds the website's interaction data unchanged
  except for the run date.
- **`ptm_pipeline_layers_raw_downloads_2026-10-18.tar.gz`** holds:
  - the two HPA 25.1 downloads;
  - Shishkova et al. 2017 Supplementary Data 1;
  - the UniProt 2026_03 FASTA as downloaded on 18 October 2026, with UniProt's response headers (which record the
    release).
- **`ptm_pipeline_RESTRICTED_cplm4_2026-10-18.tar.gz`** holds the CPLM 4.0 Methylation download and the two tables
  scripts/39 made from it. CPLM states "All Rights Reserved" and no licence, so keep this file private or deposit it
  with restricted access.
- **`PTM_Pipeline_all_branches_2026-10-20.bundle`** is a Git bundle of every branch and tag of the GitHub
  repository as of 20 October 2026: `main`, the release branch `claude/epic-franklin-oni9ft`, this branch, and the
  others. It includes the release branch's v1.1.0 tables and its own inputs (Maron 2021 supplements, ProMetheusDB,
  iPTMnet PTM master, InterPro, MobiDB, CD-CODE exports, PhosphoSitePlus terms).

## Restore

```
git clone PTM_Pipeline_all_branches_2026-10-20.bundle PTM_Pipeline   # the repository, all branches
cd PTM_Pipeline && git checkout claude/awesome-thompson-7j9n92
mkdir -p data && tar -xzf ../ptm_pipeline_interactions_raw_cache_2026-10-16.tar.gz -C data
python scripts/37b_check_staged_interactions.py                       # tables vs cache: 30 checks
python scripts/37_export_interaction_data.py --out /tmp/check.json    # rebuild the website interaction data from the cache
```

`tar` keeps each file's original date, which is the day it was downloaded.

## Upload to Zenodo

1. **Confirm the creators.** In `zenodo_public.json` and `zenodo_restricted.json`, check the creators' full names,
   and add ORCIDs (`"orcid"`) and affiliations (`"affiliation"`). The names there are copied from the paper's
   citation in the repository.
2. **Get a token.** Create a personal access token at zenodo.org (Account → Applications) with the scopes
   `deposit:write` and `deposit:actions`, then `export ZENODO_TOKEN=...`.
3. **Upload the public record:**

   ```
   python zenodo_upload.py --metadata zenodo_public.json README.md MD5SUMS SHA256SUMS \
       ptm_pipeline_inputs_tables_2026-10-20.tar.gz \
       ptm_pipeline_interactions_raw_cache_2026-10-16.tar.gz \
       ptm_pipeline_layers_raw_downloads_2026-10-18.tar.gz
   ```

   This creates a **draft**, uploads the files, checks each MD5 and sets the metadata. Review the draft on
   zenodo.org, then press Publish. Add `--sandbox` (with a sandbox.zenodo.org token) for a trial run first.
4. **Upload the restricted record (optional):**

   ```
   python zenodo_upload.py --metadata zenodo_restricted.json ptm_pipeline_RESTRICTED_cplm4_2026-10-18.tar.gz
   ```

5. **Archive the code release.** For the v1.1.0 code and release tables, copy `github_release.zenodo.json` to the
   root of the release branch as `.zenodo.json`. Turn on the repository in Zenodo's GitHub settings, then make the
   GitHub release `v1.1.0`. Zenodo archives the tagged repository and mints its DOI. Cross-link the two records
   with "is supplemented by".

The web form works too: New upload, add the files, and copy the fields from `zenodo_public.json`.

## Licences

The deposit as a whole is CC BY-NC-SA 4.0. Each source keeps its own licence, as listed in the MANIFEST files.
Checked on 20 October 2026:

| Source | Licence |
|---|---|
| UniProt, QuickGO/GO, STRING, ProtVar, AlphaFold DB, IntAct, Human Protein Atlas | CC BY 4.0 |
| Complex Portal, PDB (RCSB, ModelServer) | CC0 |
| PDBe-KB, SIFTS | EMBL-EBI terms of use (free use with attribution) |
| RBDmap (Castello et al., Mol Cell 2016) | CC BY 4.0 (open access article) |
| CARM1 knockout (Shishkova et al., Nat Commun 2017) | CC BY 4.0 |
| pepR-MS (Sha, Kuster & Trendel, bioRxiv 2025) | CC BY-NC |
| Trendel et al., Cell 2025 | CC BY-NC 4.0 (open access article) |
| CD-CODE 2.3 | CC BY-NC 4.0 (cd-code.org) |
| CPLM 4.0 | All rights reserved: restricted archive only |

## Not in this archive

- **Data for the sites v1.1.0 adds.** These are the per-site structure measurements and the PDBe-KB and AlphaFold
  DB lookups for the 873 sites and 94 proteins that release v1.1.0 adds to the published union. They have not been
  downloaded yet.
- **The rebuilt v1.1.0 release tables.** The release is to be rebuilt on UniProt 2026_03 (see
  `inputs/reference/README.md`). Its current tables are in the bundle, on branch `claude/epic-franklin-oni9ft`.
