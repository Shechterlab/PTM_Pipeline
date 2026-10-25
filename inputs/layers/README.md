# Layer inputs

These tables hold the annotation layers. The website shows them on protein and site pages; none of them adds
methylarginine sites to the union. Each layer is registered in `scripts/release_layers.py`, which lists its files (the
release build checks them against `../MANIFEST.tsv` first) and the columns, tables and tests each adds to
the release (`../../ADDING_A_SOURCE.md` says how to add one).

| Layer | Made by | Files |
|---|---|---|
| CPLM 4.0 lysine methylation | `scripts/39_stage_cplm4_lys_methylation.py` | `cplm4_human_lys_methyl_validated.tsv`, `cplm4_human_lys_methyl_rejected.tsv`, `cplm4_summary.json` |
| CARM1 knockout response (Shishkova et al. 2017) | `scripts/41_stage_carm1_ko_sites.py` | `carm1_ko_peptides.tsv`, `carm1_ko_sites.tsv`, `carm1_ko_summary.json` |
| Human Protein Atlas 25.1 subcellular locations | `scripts/42_stage_hpa_subcellular.py` | `hpa_subcellular_locations.tsv`, `hpa_summary.json` |
| Cancer cell lines and tumours: DepMap Public 24Q4 (CRISPR gene effect, changes at arginines, MTAP copy number), PRISM Repurposing 24Q2 (PRMT inhibitors) and the TCGA MC3 somatic calls at arginines (GDC open access) | `scripts/44b_extract_tcga_mc3.py` (from the MC3 public MAF, the GDC cases of TCGA and UniProt's Ensembl cross-references), then `scripts/44_stage_depmap.py`; the release build adds its columns and tests (`scripts/release_layers.py`) and `scripts/44c_export_depmap.py` writes the website's `docs/data/depmap/` from the release | `depmap_24q4_gene_effect_summary.tsv`, `depmap_24q4_arginine_changes.tsv`, `depmap_24q4_models.tsv`, `depmap_24q4_prmt_gene_effect.tsv`, `prism_24q2_prmt_inhibitors.tsv`, `tcga_mc3_arginine_changes.tsv`, `tcga_mc3_projects.tsv`, `depmap_summary.json` |
| Writer perturbations: PRMT7 knockout (Szewczyk et al. 2020), PRMT1 knockdown (Hartel et al. 2019), PRMT5 inhibitor (Lu et al. 2024), PRMT7, CARM1 and PRMT5 knockdowns and the PRMT7 inhibitor SGC3027 (Li et al. 2021) | `scripts/45_stage_prmt7_szewczyk2020.py`, `scripts/46_stage_hartel2019_prmt1_kd.py`, `scripts/51_stage_lu2024_prmt5.py`, `scripts/50_stage_li2021_perturbations.py` | `prmt7_ko_*`, `hartel2019_prmt1_kd_*`, `lu2024_prmt5_*`, `li2021_*` |
| PaxDb 6.1 protein abundance | `scripts/47_stage_paxdb.py` | `paxdb_human_integrated_abundance.tsv`, `paxdb_summary.json` |
| Reactome 97 pathways | `scripts/48_stage_reactome.py` | `reactome_human_protein_pathways.tsv`, `reactome_human_pathways.tsv`, `reactome_summary.json` |

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

- **`tcga_mc3_arginine_changes.tsv`** (274,408 rows, 9,064 tumour samples): somatic missense and stop changes at
  arginines of reviewed canonical proteins in TCGA tumours, from the MC3 public MAF of the GDC
  (`mc3.v0.2.8.PUBLIC.maf.gz`, Ellrott et al. 2018; the calls of the PanCancer Atlas). Calls are kept as MC3 kept them:
  FILTER `PASS`, and in OV and LAML also `wga` alone. The change is read on a transcript that UniProt 2026_03 assigns to
  the entry's displayed isoform, from the MAF's `all_effects` (`mapped_by` `canonical_transcript`, 95% of rows); where
  the variant falls in no such transcript, on the MAF's own transcript, placed through the UniProt entry it names or
  the gene symbol (`maf_transcript`). Either way a change counts only where the canonical sequence has an arginine at
  its position. GDC open-access data carry no licence; the GDC's data access policies (NIH Genomic Data Sharing policy)
  ask that TCGA and the MC3 publication be acknowledged and that no one try to identify participants.
- **`tcga_mc3_projects.tsv`** (33 rows): each TCGA project, with its tumour samples in the MAF (10,294 in all), those
  with a kept call, and the arginine changes extracted.

## CPLM tables are not redistributed

CPLM 4.0 states "Copyright © 2010-2026 The CUCKOO Workgroup. All Rights Reserved" and gives no licence
(cplm.biocuckoo.cn, checked 20 October 2026). For that reason:

- the two CPLM tables stay out of Git (`.gitignore`), and are marked `checksum_only` in `../MANIFEST.tsv`, which
  gives the MD5 of both;
- the CPLM summary is in Git;
- anyone who downloads `Methylation.zip` (MD5 in `raw_files.tsv`) can rebuild them with scripts/39.

The offline archive keeps them in a separate, restricted file.
