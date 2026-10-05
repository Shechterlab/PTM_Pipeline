#!/usr/bin/env Rscript
# Re-run of the 2022 enrichKEGG analysis (Maron et al. 2021 Table S5, compiled methylarginine proteins)
# with the current KEGG release, using the same call: enrichKEGG(gene, organism = "hsa",
# keyType = "ncbi-geneid", pvalueCutoff = 0.05) and KEGG's default background (all KEGG-annotated genes).
# Inputs: (a) Maron 2021 Table S5 proteins only (the 2022 input), (b) the full PTM_pipeline union.
# Gene symbols are mapped with org.Hs.eg.db instead of the 2022 GRCh37 biomaRt query.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Hs.eg.db)
})

results_root <- Sys.getenv("PTM_RESULTS_ROOT", "PTM_results")
out_dir <- file.path(results_root, "disease_gene_tests", "kegg_2022_rerun")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

union <- read.delim(file.path(results_root, "integrated", "human_arg_methyl_union_dedup_by_site.tsv"),
                    stringsAsFactors = FALSE)
inputs <- list(
  maron2021_table_s5 = unique(union$substrate_genename[grepl("Maron2021", union$source_families_exact)]),
  ptm_pipeline_union = unique(union$substrate_genename)
)

summary_rows <- list()
for (nm in names(inputs)) {
  symbols <- inputs[[nm]][inputs[[nm]] != ""]
  map <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez <- unique(map$ENTREZID)
  ek <- enrichKEGG(gene = entrez, keyType = "ncbi-geneid", organism = "hsa", pvalueCutoff = 0.05)
  ek <- setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  df <- as.data.frame(ek)
  write.table(df, file.path(out_dir, sprintf("enrichKEGG_%s.tsv", nm)), sep = "\t", quote = FALSE, row.names = FALSE)

  # Current enrichplot maps p.adjust to fill (2022 versions used colour); apply the 2022 palette to both.
  p <- dotplot(ek, showCategory = 20, label_format = 45) +
    scale_fill_gradient(low = "#E76F00", high = "#55439B") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  for (ext in c("pdf", "png", "svg")) {
    ggsave(file.path(out_dir, sprintf("enrichKEGG_%s_dotplot_2022style.%s", nm, ext)), p, width = 7.5, height = 7, dpi = 300)
  }
  summary_rows[[nm]] <- data.frame(input = nm, symbols = length(symbols), entrez_mapped = length(entrez),
                                   genes_in_kegg = sub(".*/", "", df$GeneRatio[1]), significant_pathways = nrow(df))
}
write.table(do.call(rbind, summary_rows), file.path(out_dir, "enrichKEGG_input_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
