# Inputs

Every file the releases and the website are built from, as plain tables wherever the source allows. Nothing here is
downloaded at build time: `hpc/build_release.sh` checks each file against `MANIFEST.tsv` before it runs.

| Folder | Holds |
|---|---|
| `sources/` | where the sites come from: Maron et al. 2021 Tables S2, S3 and S5, ProMetheusDB (Massignani et al. 2022), the iPTMnet-format PTM master table, and from v1.2.0 one table of localised records per study and database (`<study>_methylarginine_records.tsv`: Larsen 2016, Maron 2021, Guo 2014, Musiani 2019, Fong 2019, Li 2021, Fedoriw 2019, UniProt, Lim 2020, Wei 2020, IEDB, and the studies added in v1.2.0, Hartel 2019, Geoghegan 2015, Sylvestersen 2014, Lu 2025, Ma 2023) |
| `reference/` | the one UniProt reference: release 2026_03, reviewed human reference proteome UP000005640 (sequences, entries, accession history, release record); see its README |
| `annotation/` | InterPro 108.0 domain intervals and MobiDB disorder intervals for every reviewed entry, the PANTHER 19.0 protein classes, and the CD-CODE condensate tables |
| `layers/` | the CARM1 knockout response (Shishkova et al. 2017), Human Protein Atlas 25.1 locations, and the CPLM 4.0 summary; see its README |
| `interactions/` | the data behind the website's interaction views (STRING, Complex Portal, GO, PDBe-KB and SIFTS, PDB contacts, ProtVar, AlphaFold DB, IntAct, RNA and DNA crosslinks, RBDmap, CD-CODE 2.3), one table per source; see its README |
| `disease/` | DisGeNET term-gene sets as bundled with DOSE 4.4.0, org.Hs.eg.db 3.22.0's gene symbol to Entrez map, and the Human Disease Ontology (release 2026-09-30) with the gene sets DOSE enrichDO tests, for the disease enrichment |
| `licences/` | the PhosphoSitePlus terms as of July 2021 (web archive capture), when the PhosphoSitePlus records in Maron et al. 2021 Table S5 were downloaded |

## How each folder was staged

Staging is done once, before a build; the build itself reads only these files. Every file's own row in
`MANIFEST.tsv` says exactly how it was made (`how_made_or_checked`) and from what (`derived_from`).

| Folder | Staged with |
|---|---|
| `sources/` | the journal supplements as published (checked against the publishers' files), and the lab's PTM master table. The study tables: `scripts/00_stage_study_sites.py` reads each study's own tables (supplementary tables; for Maron 2021, Tables S2 and S3 with the lab's Proteome Discoverer exports of the same immunoprecipitations; URLs, MD5s and the command in `MANIFEST.tsv`) and writes one record per reported peptide and arginine whose methyl group the study localises, with the filters of the study's methods. Raw downloads go to the Zenodo archive, except files whose terms do not allow it, of which only derived records are redistributed (`checksum_only`): the Elsevier supplements of Fong 2019 and Fedoriw 2019, the Science Signaling tables of Larsen 2016 and Musiani 2019 (AAAS), and the unpublished Maron 2021 Proteome Discoverer exports |
| `reference/` | UniProt REST downloads (query and release in `uniprot_release.json`); the accession history with `scripts/00_stage_context_sources.py` |
| `annotation/` | InterPro: `scripts/00_stage_context_sources.py --interpro-source bulk` on InterPro 108.0 `protein2ipr.dat.gz`. MobiDB: `scripts/00_stage_mobidb_disorder.py`. Both moved to the 2026_03 sequences by `scripts/00_lift_intervals_to_reference.py`. PANTHER: `scripts/00_stage_panther_classes.py` on the three PANTHER 19.0 downloads (kept at their release URLs). UniProt methylarginine annotations: `scripts/00_stage_uniprot_methylarginine.py` (UniProt REST, release 2026_03). CD-CODE: exports as downloaded |
| `layers/` | downloads by `scripts/00_fetch_static_sources.py`, then `scripts/41_stage_carm1_ko_sites.py`, `scripts/42_stage_hpa_subcellular.py` and `scripts/39_stage_cplm4_lys_methylation.py` |
| `interactions/` | `scripts/37_export_interaction_data.py` downloads into a cache, and `scripts/37a_stage_interactions.py` turns the cache into these tables. For proteins and sites a release adds, `scripts/37 --tables` fetches what the tables lack into a top-up cache, and `scripts/37a_stage_interactions.py --topup` stages it as dated tables next to the others |
| `disease/` | `scripts/00_export_disease_tables.R` (Bioconductor DOSE 4.4.0 and org.Hs.eg.db 3.22.0, from bioconda) and `scripts/00_export_disease_ontology.R` (DOSE's Disease Ontology data, HDO.sqlite from YuLab) |

## MANIFEST.tsv

One row per file: path, size, MD5, rows, `redistribution`, content, source URL, version or date, citation, licence,
`derived_from` (the raw file a table was made from) and how it was made or checked.

`redistribution` says where a file is:

- `git`: in this repository;
- `zenodo`: in the Zenodo archive, not in Git (the raw downloads, and the inputs only release v1.0.0 used);
- `checksum_only`: not redistributed, because the licence does not allow it (CPLM 4.0, "All Rights Reserved") or
  because the file stays at a permanent public URL (InterPro 108.0's 17 GB `protein2ipr.dat.gz`); the MD5 lets
  anyone check their own copy.

`interactions/raw_cache_files.tsv` and `layers/raw_files.tsv` give the MD5 of every raw download behind those
tables. Independent re-downloads of the bulk files on 20 October 2026 had the same MD5 as the copies the website was
built from.

## Notes

- **Names** carry the source version or the download date.
- **InterPro and MobiDB coordinates** were computed on older UniProt sequences. 13 entries changed sequence in
  2026_03. Their intervals are moved to 2026_03 coordinates where the change is an N-terminal extension, and are kept
  where only single residues changed; the six other entries lose their intervals (see the annotation rows of
  `MANIFEST.tsv`).
- **Licences** are those of each source (`MANIFEST.tsv`). PhosphoSitePlus records inside Maron Table S5 are CC
  BY-NC-SA 3.0 Unported, per `licences/`.
