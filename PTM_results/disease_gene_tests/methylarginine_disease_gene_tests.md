# Methylarginine in ALS genes: disease-ontology, ClinGen and TDP-43 tests

Scripts: `scripts/36_build_disease_test_gene_table.py` → `37_methylarg_dose_kegg_enrichment.R` → `38_methylarg_clingen_tdp43_tests.py` → `39_plot_methylarg_disease_figures.py`.
Run date 2026-10-05. Versions: R 4.6.1, DOSE 4.6.0, clusterProfiler 4.20.0 (`dose_versions.txt`).

## Inputs

| Input | Source | Count |
|---|---|---|
| Methylarginine genes (any evidence) | `PTM_results/integrated/human_arg_methyl_union_dedup_by_site.tsv` (Maron 2021 Table S5, PRoMetheusDB, curated PTM master). 4,613 proteins mapped by reviewed accession, 431 by gene symbol (TrEMBL accessions), 46 unmapped | **4,687 genes** |
| Gene universe | UniProt reviewed human with Entrez ID (downloaded 2026-10-05) | 19,113 |
| MS-detectable background | PaxDb 5.0 cell-line proteomes (Geiger 2012 iBAQ and spectral counts, 11 lines; Nagaraj 2011 ×3; Kulak 2014; Müller 2020), 11,320 genes, plus every methylarginine gene | 11,945 |
| Covariates | RNA binding = GO:0003723 incl. children (org.Hs.eg.db); log10 protein length; log10 PaxDb whole-organism abundance | — |
| Disease genes | ClinGen gene–disease validity (downloaded 2026-10-05), Definitive/Strong/Moderate; 42 curation panels with ≥10 genes | 2,471 genes |
| TDP-43 targets | `Shechterlab/cryptic_exon` `cryptic_sf3b2_screen/results/external_tdp43_compendium_gene_summary.hg38.tsv` | 3,385 genes |

PXD076372 was deliberately not used. "Heavily methylated" = ≥3 sites in the union.

## Results

### 1. DOSE Disease Ontology (MS-detectable background) — `fig_disease_ontology_bar`

| Term | Rank (of 864) | Genes | Fold | p | FDR |
|---|---|---|---|---|---|
| Motor neuron disease | 8 | 78 of 137 | 1.33 | 0.0006 | 0.06 |
| ALS | 18 | 55 of 97 | 1.32 | 0.004 | 0.19 |

Top terms are autosomal dominant intellectual developmental disorder, intellectual disability and specific developmental disorder. Against the full reviewed proteome, motor neuron disease ranks 6 of 1,045 (FDR 0.0008) and ALS 20 (FDR 0.014). In DisGeNET (MS-detectable), Motor Neuron Disease ranks 35 of 7,085 (FDR 0.01), and the top terms are neurodevelopmental/dysmorphic. Frontotemporal dementia is not enriched in either ontology with this gene set.

### 2. KEGG (MS-detectable background) — `fig_kegg_simple_dotplot`

Spliceosome is first (113 genes, FDR 6e-18), followed by chromatin remodeling, ribosome, mRNA surveillance and nucleocytoplasmic transport. The KEGG ALS pathway ranks 11 of 344 (165 genes, fold 1.26, FDR 0.001). This is a simplified version of the 2022 enrichKEGG dot plot: one colour, fold enrichment on x, and gene count and FDR given as text.

### 3. ClinGen panels — `clingen_panel_methylarginine.tsv`

- ALS Spectrum Disorders: **21 of 29 genes (72%) methylated, rank 1 of 42 panels**. For comparison, 35% of all curated disease genes and 25% of the reviewed proteome are methylated.
- Adjusted OR 3.04 (95% CI 1.25–7.39), p = 0.014.
- Methylated ALS genes: ALS2, ANXA11, AR, ATXN2, CHMP2B, DCTN1, FUS, HNRNPA2B1, KIF5A, LRP12, MATR3, NEK1, PFN1, SOD1, SPTLC2, SQSTM1, TARDBP, TUBA4A, UBQLN2, VAPB, VCP.
- Heavily methylated (≥3 sites): ALS ranks 5 by fraction (6 of 29), but the adjusted OR is 0.72 (p = 0.54), so it is **not enriched**.
- With 42 panels, p = 0.014 does not survive correction. Report the result as a rank (1 of 42), not as a corrected p.

### 4. TDP-43 compendium — `fig_als_vs_tdp43_forest`, `tdp43_compendium_methylarginine.tsv`

Adjusted ORs on the MS-detectable background:

| Gene set | Methylated | Adjusted OR (95% CI) | p |
|---|---|---|---|
| Inherited ALS genes (ClinGen), comparator | 21/29 | 2.52 (1.06–5.96) | 0.036 |
| All compendium host genes | 908/2,348 | 0.89 (0.80–0.98) | 0.017 (slight depletion) |
| TDP-43-sensitive cryptic events | 155/358 | 1.04 (0.83–1.30) | 0.72 |
| NMD-rescued / dark cryptic events | 272/681 | 0.90 (0.76–1.06) | 0.21 |
| Down after TDP-43 knockdown (≥2 datasets) | 326/840 | 0.94 (0.81–1.09) | 0.41 |
| Pathology-discriminating events (AUC ≥ 0.8) | 18/29 | 2.18 (1.01–4.72) | 0.048 |

- The large TDP-43 target sets are not enriched.
- The small AUC ≥ 0.8 set is borderline: 18 of 29 genes, mostly with 1–3 sites (AARS1, DOCK1, ARHGAP32, PFKP, IK, …). The canonical neuronal targets STMN2, UNC13A, KCNQ2, CAMK2B and CACNA1E carry no sites.
- This is one of five TDP-43 sets tested and is not corrected for multiplicity. Treat it as unconfirmed.

## Interpretation

Methylarginine marks the inherited ALS genes, which are largely RNA-binding and proteostasis proteins (FUS, TAF15, EWSR1, HNRNPA1/A2B1, MATR3, VCP, …). This holds beyond RNA binding, length and abundance. It does not mark the bulk of the RNAs that TDP-43 mis-splices downstream. The claim is genetic and made at gene-set level; it does not imply that methylarginine causes disease.

## Differences from the earlier note (MS-verified PXD076372 set, DOSE 4.4)

- **Motor neuron disease rank.** It ranked 1st with the 931 MS-verified proteins. With the 4,687 any-evidence genes it ranks 8th, behind developmental terms; the earlier note already said the any-evidence set ranks lower.
- **Background size.** The MS-detectable background is 11,945 here, against about 11,000 earlier.
- **ClinGen ALS result.** It reproduces: 21/29, rank 1, OR 3.04 here against 3.2 earlier.
- **Heavy-methylation rank.** ALS ranks 5 here, against 19 earlier; neither result is enriched after adjustment.
- **TDP-43 gene sets.** These come from the external compendium, not from Ma 2022 or Brown 2022 directly.

## ALS talk slides (`als_talk_figures/`, script 41) and KEGG 2022 re-run (`kegg_2022_rerun/`, script 40)

Counts are proteins: the union's 5,090 methylated proteins map to 4,687 genes. Ontology tools annotate genes, so each term count is the number of annotated proteins that carry a methylarginine site.

| Slide | Message | Numbers |
|---|---|---|
| `als_clingen_slide` | 72% of ALS-causing proteins are arginine-methylated | 21/29 ALS vs 35% of disease proteins vs 25% of all human proteins; rank 1 of 42 panels; adjusted OR 3.0 |
| `als_kegg_slide` | ALS is the #2 KEGG pathway | 164 of 369 pathway proteins methylated, 1.8-fold, FDR 6e-15 (KEGG default background) |
| `als_do_slide` | Motor neuron disease ranks #6 of 1,045 DO terms | 78 of 167, FDR 8e-4 (reviewed proteome background) |
| `als_protein_map_slide` | Sites sit in RGG / Gly-rich regions | FUS 30, EWSR1 35, TAF15 41, TDP-43 2 (R275, R293), hnRNPA1 20, hnRNPA2B1 15, MATR3 9 (canonical positions only) |

**KEGG 2022 re-run.** This is the same `enrichKEGG(organism = "hsa", pvalueCutoff = 0.05)` call with the default KEGG background, run on the current KEGG release.

| Input | Pathways of neurodegeneration | ALS |
|---|---|---|
| Maron 2021 Table S5 alone | #1 by GeneRatio, 173 proteins | #2 by GeneRatio, 162 proteins, FDR 7e-15 |
| Full union | #1 by GeneRatio, 176 proteins | #2 by GeneRatio, 164 proteins, FDR 6e-15; #2 by p value |

**Caveats for the speaker.**
- The KEGG and DO slides use whole-genome backgrounds, as the 2022 analysis did. Restricting to MS-detectable proteins weakens them:
  - KEGG ALS: rank 11, FDR 0.001.
  - DO motor neuron disease: rank 8, FDR 0.06.
- The KEGG ALS map (hsa05014) is large, at 369 proteins. It includes proteasome, mitochondrial and cytoskeletal proteins, not only ALS genes.
- The ClinGen result is the least background-sensitive: rank 1 of 42 panels, with covariate adjustment.
- Site counts on the protein map use canonical UniProt positions only. Rows mapped to TrEMBL or isoform accessions are excluded, so EWSR1, hnRNPA1 and MATR3 show fewer sites here than the union holds under the gene name.
