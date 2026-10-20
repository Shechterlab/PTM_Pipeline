# Frozen inputs

Every file the methylarginine releases and the website are built from, exactly as downloaded (three large IntAct tables are stored gzip-compressed; the manifest gives the original's checksum). `MANIFEST.tsv` gives each file's size, MD5,
former name, content, source, version, citation, licence and how it was checked; `hpc/build_release.sh` refuses to
run when a file differs from it.

| Folder | Holds |
|---|---|
| `sources/` | the site sources: Maron et al. 2021 Tables S2, S3 and S5, ProMetheusDB (Massignani et al. 2022) and the iPTMnet-format PTM master table |
| `reference/` | the UniProtKB reviewed human sequences every position refers to (staged 2026-04-07), their entry table, and the UniProt 2026_03 accession history |
| `annotation/` | per-protein and per-site context: MobiDB disorder intervals, InterPro 108.0 domain intervals for all reviewed entries, CD-CODE condensates, CPLM 4.0 lysine methylation, the CARM1 knockout response (Shishkova et al. 2017) and Human Protein Atlas locations |
| `interactions/` | what the website's interaction view is built from: STRING, Complex Portal, SIFTS, ProtVar, IntAct features, RNA and DNA crosslinks (pepR-MS, Trendel et al. 2025) and RBDmap |
| `licences/` | the PhosphoSitePlus terms as of July 2021 (web archive capture), when the PhosphoSitePlus records in Maron et al. 2021 Table S5 were downloaded |

The dates in the names are the staging or export dates; UniProt did not record a release number for the 2026-04-07
sequences (they predate 2026_03; see `DATA_LINEAGE_REVIEW.md`). Each file keeps its own licence (`MANIFEST.tsv`);
PhosphoSitePlus records inside Maron Table S5 are CC BY-NC-SA 3.0 Unported, per `licences/`.
