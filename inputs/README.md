# Inputs

Every file the releases and the website are built from, as plain tables wherever the source allows. Nothing here is
downloaded at build time: `hpc/build_release.sh` checks each file against `MANIFEST.tsv` before it runs.

| Folder | Holds |
|---|---|
| `sources/` | where the sites come from: Maron et al. 2021 Tables S2, S3 and S5, ProMetheusDB (Massignani et al. 2022) and the iPTMnet-format PTM master table |
| `reference/` | the one UniProt reference: release 2026_03, reviewed human reference proteome UP000005640 (sequences, entries, accession history, release record); see its README |
| `annotation/` | InterPro 108.0 domain intervals and MobiDB disorder intervals for every reviewed entry, and the CD-CODE condensate tables |
| `layers/` | the CARM1 knockout response (Shishkova et al. 2017), Human Protein Atlas 25.1 locations, and the CPLM 4.0 summary; see its README |
| `interactions/` | the data behind the website's interaction views (STRING, Complex Portal, GO, PDBe-KB and SIFTS, PDB contacts, ProtVar, AlphaFold DB, IntAct, RNA and DNA crosslinks, RBDmap, CD-CODE 2.3), one table per source; see its README |
| `licences/` | the PhosphoSitePlus terms as of July 2021 (web archive capture), when the PhosphoSitePlus records in Maron et al. 2021 Table S5 were downloaded |

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
