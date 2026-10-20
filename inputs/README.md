# Frozen inputs

Every file the methylarginine release is built from, exactly as used. `MANIFEST.tsv` gives each file's size, MD5,
former name, content, source, version, citation, licence and how it was checked; `hpc/build_release.sh` refuses to
run when a file differs from it.

| Folder | Holds |
|---|---|
| `sources/` | the site sources: Maron et al. 2021 Tables S2, S3 and S5, ProMetheusDB (Massignani et al. 2022) and the iPTMnet-format PTM master table |
| `reference/` | the UniProtKB reviewed human sequences every position refers to (staged 2026-04-07), their entry table, and the UniProt 2026_03 accession history |
| `annotation/` | MobiDB disorder intervals, InterPro domain intervals and CD-CODE condensate tables |

The dates in the names are the staging or export dates; UniProt did not record a release number for the 2026-04-07
sequences (they predate 2026_03; see `DATA_LINEAGE_REVIEW.md`). Each file keeps its own licence (`MANIFEST.tsv`);
PhosphoSitePlus records inside Maron Table S5 are CC BY-NC-SA 3.0.
