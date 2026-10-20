# Provenance

What the methylarginine site union is made from, and where each version is.

## Versions

| Version | Sites | Where |
|---|---|---|
| v1.0.0 | 13,719 sites on 5,090 accessions | Published with DeAngelo, Silverstein & Shechter, Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011. Tables in `releases/v1.0.0/`; the full analysis output at git tag `v1.0.0` and in branch `v1.0.0_archive`; the website as published in `docs/v1.0.0/`. |
| v1.1.0 | 13,358 sites on 4,722 proteins, and 10 pending sites | `releases/v1.1.0/`, and the website in `docs/`. Same sources and methods, mapped to UniProt release 2026_03, with corrections to how some source records were placed on UniProt sequences (`releases/v1.1.0/CHANGES.md`). |

## Sources of the sites

28,960 source records from three tables (`inputs/sources/`):

| Table | Records | What it is |
|---|---:|---|
| Maron et al. 2021, Table S5 | 26,515 | Methylarginine sites compiled from ten studies, PhosphoSitePlus among them (downloaded 22 May 2021) |
| ProMetheusDB (Massignani et al. 2022), mmc2 | 1,598 | Methylarginine sites with methylation states |
| iPTMnet-format PTM master table | 847 | Arginine methylation rows from UniProt (762) and IEDB (85) |

Maron et al. 2021 Tables S2 and S3 (immunoprecipitation and PRMT inhibitors) and the CARM1 knockout data of Shishkova
et al. 2017 annotate sites but add none.

## Reference and annotation

- **Sequences:** UniProt release 2026_03, the 20,416 reviewed entries of the human reference proteome UP000005640
  (`inputs/reference/`).
- **Context:** InterPro 108.0 domains and repeats, MobiDB predicted disorder (export of 2026-03-30), CD-CODE
  condensates, Human Protein Atlas 25.1 locations (`inputs/annotation/`, `inputs/layers/`).
- **Interactions:** STRING 12.5, Complex Portal, QuickGO, PDBe-KB and SIFTS, PDB contacts, ProtVar, AlphaFold DB
  complexes, IntAct, pepR-MS, RBDmap, Trendel et al. 2025 (`inputs/interactions/`).

`inputs/MANIFEST.tsv` lists every input file with its source, version or download date, licence and MD5, and
`inputs/README.md` says how each folder was staged.

## How records become sites

Each record is placed on a reviewed UniProt sequence from its own evidence: a window centred on the site, a peptide,
or a position (`scripts/02_remap_arg_methyl_to_canonical.py`). The build stops if any record is placed against its
evidence (`scripts/qc_remap_evidence_check.py`) or any site is not an arginine. `bash hpc/build_release.sh VERSION`
runs it all from `inputs/`.

## Raw downloads

The raw files behind the tables in `inputs/` (journal supplements, bulk downloads, API answers) are in the Zenodo
archive of the release; their MD5s are in `inputs/MANIFEST.tsv`, `inputs/layers/raw_files.tsv` and
`inputs/interactions/raw_cache_files.tsv`. CPLM 4.0 (all rights reserved) is not redistributed.
