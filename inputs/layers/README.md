# Layer inputs

These tables hold four annotation layers. The website shows them on protein and site pages; none of them adds
methylarginine sites to the union.

| Layer | Made by | Files |
|---|---|---|
| CPLM 4.0 lysine methylation | `scripts/39_stage_cplm4_lys_methylation.py` | `cplm4_human_lys_methyl_validated.tsv`, `cplm4_human_lys_methyl_rejected.tsv`, `cplm4_summary.json` |
| CARM1 knockout response (Shishkova et al. 2017) | `scripts/41_stage_carm1_ko_sites.py` | `carm1_ko_peptides.tsv`, `carm1_ko_sites.tsv`, `carm1_ko_summary.json` |
| Human Protein Atlas 25.1 subcellular locations | `scripts/42_stage_hpa_subcellular.py` | `hpa_subcellular_locations.tsv`, `hpa_summary.json` |
| Cancer cell lines and tumours: DepMap Public 24Q4 (CRISPR gene effect, changes at arginines, MTAP copy number), PRISM Repurposing 24Q2 (PRMT inhibitors) and TCGA PanCancer Atlas changes at arginines (cBioPortal) | `scripts/44b_fetch_tcga_arginine_changes.py` (network), then `scripts/44_stage_depmap.py`; `scripts/44c_export_depmap.py` writes the website's `docs/data/depmap/` | `depmap_24q4_gene_effect_summary.tsv`, `depmap_24q4_arginine_changes.tsv`, `depmap_24q4_models.tsv`, `depmap_24q4_prmt_gene_effect.tsv`, `prism_24q2_prmt_inhibitors.tsv`, `tcga_pancan_arginine_changes.tsv`, `tcga_pancan_studies.tsv`, `depmap_summary.json` |

- `../MANIFEST.tsv` (one manifest for all of `inputs/`) describes each file: rows, MD5, what it holds, source, version, citation, licence and how it
  was made.
- `raw_files.tsv` lists the downloads the tables were made from, with URL, size, MD5, SHA-256 and download date.
- The summary JSON files are the scripts' own summaries, with paths shortened to file names. Their input
  checksums equal those stored in the website data (`docs/data/methylarg_sites.json`, `meta.cplm4`, `meta.hpa`,
  `meta.carm1_ko`).

## Checks

On 20 October 2026 the three scripts were run again from the downloads in `raw_files.tsv`. All five tables came
out byte-identical, both times:

- with the FASTA the scripts originally used, the UniProt 2026_03 reviewed human reference proteome
  (20,416 sequences). That file is now in Git as `../reference/uniprot_human_reviewed_2026_03.fasta`;
- with `../interactions/uniprot_human_reviewed_2026-10-16.fasta`, also in Git.

On 21 October 2026 the HPA tables were made again after a fix to scripts/42. Where HPA maps one accession to
several genes, the script had kept the reliability and main/additional flag of the first gene listed; it now keeps
the best reliability and counts a location as main when any of the genes has it as main. Four rows changed.

## What each table holds

- **`cplm4_human_lys_methyl_validated.tsv`** (5,356 sites on 2,939 proteins): CPLM 4.0 human lysine methylation
  sites whose CPLM sequence matches the canonical UniProt sequence. The match is either an identical sequence, or
  the same ±7-residue window around the lysine. Nothing is remapped. The CPLM id, evidence label (`Exp.`, `Dat.`,
  `Exp./Dat.`) and references are given as CPLM gives them.
- **`cplm4_human_lys_methyl_rejected.tsv`** (559 rows): the rows that failed, with the reason:
  - 494: the accession is not in the reference proteome;
  - 60: no lysine at that position in the canonical sequence;
  - 5: the window differs.
- **`carm1_ko_peptides.tsv`** (515 peptides): ADMA peptides quantified in parental and CARM1 knockout MCF7 and
  MDA-MB-231 cells. Each peptide is placed on the canonical sequence and carries its log2 KO/WT, p value and call:
  - down: log2 ≤ −1 and p < 0.05;
  - up: log2 ≥ 1 and p < 0.05.
- **`carm1_ko_sites.tsv`** (603 rows, 431 sites): one row per methylarginine and cell line, with:
  - median log2 KO/WT;
  - lowest p;
  - peptide counts;
  - the call (down, up, unchanged, mixed).
- **`hpa_subcellular_locations.tsv`** (27,150 rows on 13,301 proteins): one row per protein and HPA location,
  with:
  - whether it is a main (1) or additional (0) location;
  - the HPA reliability (Enhanced, Supported, Approved, Uncertain).

## CPLM tables are not redistributed

CPLM 4.0 states "Copyright © 2010-2026 The CUCKOO Workgroup. All Rights Reserved" and gives no licence
(cplm.biocuckoo.cn, checked 20 October 2026). For that reason:

- the two CPLM tables stay out of Git (`.gitignore`), and are marked `checksum_only` in `../MANIFEST.tsv`, which
  gives the MD5 of both;
- the CPLM summary is in Git;
- anyone who downloads `Methylation.zip` (MD5 in `raw_files.tsv`) can rebuild them with scripts/39.

The offline archive keeps them in a separate, restricted file.
