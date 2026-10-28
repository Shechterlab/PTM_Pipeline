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

The card's RNA and DNA layers (`pilot_layers.py`, written to `docs/data/subcellular_na.json`) place the transcripts and
the genes of proteins, each gene on its reviewed UniProt entry by gene symbol, and compare the genes of methylated
proteins with the other genes each dataset places, within deciles of their expression in the same cells.

| Layer | Map | File read | Source |
|---|---|---|---|
| RNA | APEX-seq (Fazal et al. 2019, HEK293T, poly(A)) | `GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz`, `..._counts_all_data_polyA_hek293t.txt.gz` | GEO GSE116008 |
| RNA | Nuclear-body APEX2-seq (Barutcu et al. 2022, Flp-In HEK293, ribosome-depleted) | `GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz` | GEO GSE176439 |
| DNA | TSA-seq (Chen et al. 2018, K562): SON and lamin B, two replicates each | `GSE81553_SON_TSA-Seq_Condition2*.bw`, `GSE81553_LaminB_TSA-Seq_Conditoin2*.bw`; GENCODE v19 genes | GEO GSE81553; EBI GENCODE |
| DNA | expression for the deciles: K562 total RNA-seq | ENCODE `ENCFF611MXW.tsv`, `ENCFF384BFE.tsv` (ENCSR109IQO) | ENCODE |

Rules: APEX-seq, enriched over the compartment's control at FDR < 0.05 (the study's DESeq2 table; ER lumen has no
enriched transcript and is left out); nuclear bodies, labelled at least twice as much as by nucleoplasmic APEX2 in both
replicates (speckles two of three baits, PML bodies both PML and SP100; the nucleolar FBL bait places only 3 transcripts
and is not shown); TSA-seq, signal over the gene (at least 20 kb) in the top tenth of the genome's 20 kb bins.
RD-SPRITE has no human data (GEO GSE151515 is mouse ES cells); human SPRITE (GM12878, GSE114242) is DNA only.

The inputs are not in the repository yet; v1.5.0 would stage them into `inputs/` with checksums, as for every layer.

    python pilot_hcm.py REPO HCM_DOWNLOADS OUT
    python pilot_oip.py REPO ORGANELLE_IP_REPO OUT
    python pilot_chrom.py REPO embj0033-0648-sd8.xlsx OUT
    python pilot_layers.py REPO DOWNLOADS OUT       # needs pyBigWig
    python build_preview.py REPO OUT docs/data/subcellular.json   # also writes docs/data/subcellular_na.json
