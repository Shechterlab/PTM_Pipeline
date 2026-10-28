# Where in the cell: pilot (preview of v1.5.0, not part of a release)

Each map is tested against the proteins it places: methylated proteins (release v1.4.1) per compartment, one-sided
hypergeometric with fold enrichment, and the Mantel-Haenszel odds ratio within PaxDb abundance deciles of the same
proteins (`abundance_strata.py`). `build_preview.py` writes `docs/data/subcellular.json`, which the website's "Where in
the cell" card (preview) and the protein page's "Where the maps place it" read. The file carries each protein's PaxDb
decile among the proteins each map places, and the HPA locations of every HPA-profiled protein, so that the card repeats
the tests within the deciles for a selection (Mantel-Haenszel in the browser; with all methylated proteins as the
selection it gives the pilots' counts, odds ratios and confidence intervals).

| Map | File read | Source and licence |
|---|---|---|
| Human Cell Map (Go et al. 2021, HEK293) | `preys-v1.txt`, `not-significant-v1.txt` | https://humancellmap.org/resources/downloads (CC BY 4.0) |
| Organelle IP (Hein et al. 2025, HEK293T) | `data/labels/cluster_annotation_Dec6.csv` | github.com/czbiohub-sf/Organelle_IP_analyses_and_figures (BSD-3) |
| Chromatin probability (Kustatscher et al. 2014) | `embj0033-0648-sd8.xlsx` (Supplementary Table 1) | EMBO J 33:648, open access (CC BY) |
| Human Protein Atlas 25.1 | release statistics v1.4.1 (`hpa_location`) | already in the release (CC BY 4.0) |
| Speckle, lamina and PML-body proteomes (Barutcu et al. 2022, HEK293) | Mol Cell Table S1 (`mmc2.xlsx`) | Elsevier supplementary table; calls derived from it |

A fifth protein map, the speckle, lamina and PML-body proximity proteomes of Barutcu et al. 2022 (APEX2-MS, Flp-In
HEK293; Mol Cell Table S1, SAINT; raw data MassIVE MSV000087649), comes from `pilot_apexms.py`: a protein belongs where
SAINT's average probability is at least 0.9 (speckles: with two of the three markers SRSF1, SRSF7, RNPS1; lamina: LMNA
or LMNB1; PML bodies: SP100); the proteins it places are every prey detected.

The card's RNA and DNA layers (`pilot_layers.py`, written to `docs/data/subcellular_na.json`) ask whether methylated
proteins bind RNAs that gather in particular compartments more often than other RNA-binding proteins do. Each protein
with eCLIP data (ENCODE and ENCORE, from `inputs/nucleic_acid/`: the genes with most reproducible peaks in each
experiment, pooled over a protein's experiments) is a unit. Its bound RNAs gather in a compartment when an RNA map
(APEX-seq, nuclear-body APEX2-seq) places them there more often than the other RNAs it places, within deciles of
expression in the map's cells crossed with deciles of gene length (GENCODE v19; Mantel-Haenszel; BH over every protein
and compartment of the map; q <= 0.05, odds ratio above 1); the DNA maps (TSA-seq, Hi-C subcompartments) place the genes of those RNAs the same way. A map tests the
proteins with at least 20 bound RNAs it places, and compartments where it places fewer than 20 RNAs (genes) are not
tested. The methylated proteins are then compared with the other eCLIP'd proteins a map tests, within PaxDb abundance
deciles (`lib.tests`, as for the protein maps; compartments where fewer than five proteins' RNAs gather are not
compared). A call depends on how many bound RNAs a map places, so each compartment also compares the size of every
protein's odds ratio, called or not: methylated against the others, two-sided Mann-Whitney U, BH per map. Genes are
joined by Ensembl gene ID (nuclear-body counts, by a symbol naming a single gene).

Gene length has to be matched as well as expression: eCLIP target lists favour long genes, whose introns collect more
peaks, and compartments differ in it. The RNAs nuclear-body APEX2-seq places at the lamina come from genes of median
319 kb (71% over 200 kb, against 9% of all it places), its speckle RNAs from genes of 12 kb. Matched on expression
alone, 24 proteins' bound RNAs "gathered" at the lamina; matched on both, one.

Of the 230 eCLIP'd proteins, 186 carry a methylarginine; the maps test 219-220 (178-179 methylated), so every comparison
rests on about 41 proteins without one. No compartment separates the two groups (lowest q 0.11: APEX-seq lamina, odds
ratio 2.6), and neither do the sizes of their odds ratios.

Exploratory, not on the page: within the methylated eCLIP'd proteins, those with a methylarginine at a residue-level RNA
contact (an RNA crosslink at the residue, an RBR-ID residue, a CAPRI crosslinked peptide or a PDB contact; 65 of 179)
bind RNAs from speckle-associated genes less often than the rest (TSA-seq odds ratio 0.38, q 0.01; Hi-C A1 0.34,
q 0.007, within PaxDb deciles). Against proteins whose residue-level RNA contacts are all at other residues (82), the
difference shrinks (0.57, q 0.19; 0.50, q 0.12): it belongs mostly to well-characterised sequence-specific RNA-binding
proteins, not to where their methylarginines sit. Proteins with a methylarginine anywhere in RNA-binding evidence
(including RBDmap fragments and pCLAP and CAPRI peptides; 114 of 179) carry a median of 7 methylarginines against 1
for the rest, so that split mostly sorts proteins by how many sites they have.

| Layer | Map | File read | Source |
|---|---|---|---|
| RNA | APEX-seq (Fazal et al. 2019, HEK293T, poly(A)) | `GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz`, `..._counts_all_data_polyA_hek293t.txt.gz` | GEO GSE116008 |
| RNA | Nuclear-body APEX2-seq (Barutcu et al. 2022, Flp-In HEK293, ribosome-depleted) | `GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz` | GEO GSE176439 |
| DNA | TSA-seq (Chen et al. 2018, K562): SON and lamin B, two replicates each | `GSE81553_SON_TSA-Seq_Condition2*.bw`, `GSE81553_LaminB_TSA-Seq_Conditoin2*.bw`; GENCODE v19 genes | GEO GSE81553; EBI GENCODE |
| DNA | genes placed and their expression: K562 total RNA-seq (TPM >= 1) | ENCODE `ENCFF611MXW.tsv`, `ENCFF384BFE.tsv` (ENCSR109IQO) | ENCODE |
| DNA | Hi-C subcompartments (Rao et al. 2014, GM12878): A1, A2, B1, B2, B3, B4 | `GSE63525_GM12878_subcompartments.bed.gz` | GEO GSE63525 |
| DNA | genes placed and their expression: GM12878 total RNA-seq (TPM >= 1) | ENCODE `ENCFF723ICA.tsv`, `ENCFF362RMV.tsv` (ENCSR000AEE) | ENCODE |
| both | RNAs bound: eCLIP targets | `encode_eclip_target_genes.tsv`, `encore_eclip_target_genes.tsv` | already in `inputs/` |

Rules: APEX-seq, enriched over the compartment's control at FDR < 0.05 (the study's DESeq2 table; ER lumen has no
enriched RNA and is left out); nuclear bodies, labelled at least twice as much as by nucleoplasmic APEX2 in both
replicates (speckles two of three markers, PML bodies both PML and SP100); TSA-seq, signal over the gene (at least 20 kb)
in the top tenth of the genome's 20 kb bins; Hi-C, the subcompartment at the middle of the gene (GENCODE v19, hg19, as
the bed file).
RD-SPRITE has no human data (GEO GSE151515 is mouse ES cells); human SPRITE (GM12878, GSE114242) is DNA only.

The inputs are not in the repository yet; v1.5.0 would stage them into `inputs/` with checksums, as for every layer.

    python pilot_hcm.py REPO HCM_DOWNLOADS OUT
    python pilot_oip.py REPO ORGANELLE_IP_REPO OUT
    python pilot_chrom.py REPO embj0033-0648-sd8.xlsx OUT
    python pilot_apexms.py REPO mmc2.xlsx OUT       # Barutcu et al. 2022 Table S1
    python pilot_layers.py REPO DOWNLOADS OUT       # needs pyBigWig
    python build_preview.py REPO OUT docs/data/subcellular.json   # also writes docs/data/subcellular_na.json
