#!/usr/bin/env Rscript
# Disease Ontology (enrichDO), DisGeNET (enrichDGN) and KEGG (enrichKEGG) over-representation
# of methylarginine genes, against two backgrounds:
#   ms_detectable : genes detected in PaxDb 5.0 cell-line proteomes, plus all methylarginine genes
#   reviewed      : all reviewed human UniProt entries with an Entrez ID
# Also exports GO:0003723 (RNA binding, including child terms) membership for the covariate models.

suppressPackageStartupMessages({
  library(DOSE)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
results_root <- Sys.getenv("PTM_RESULTS_ROOT", "PTM_results")
out_dir <- if (length(args) >= 1) args[[1]] else file.path(results_root, "disease_gene_tests")
gene_table <- read.delim(file.path(out_dir, "gene_table.tsv"), stringsAsFactors = FALSE)
gene_table$entrez <- as.character(gene_table$entrez)

query <- gene_table$entrez[gene_table$methyl_any == "True"]
backgrounds <- list(
  ms_detectable = gene_table$entrez[gene_table$ms_detectable == "True"],
  reviewed = gene_table$entrez
)

rna <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0003723", keytype = "GOALL", columns = "ENTREZID")
gene_table$rna_binding <- gene_table$entrez %in% unique(rna$ENTREZID)
write.table(gene_table[, c("entrez", "rna_binding")], file.path(out_dir, "gene_rna_binding_go0003723.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

run_one <- function(fun, bg) {
  res <- fun(gene = query, universe = bg, pvalueCutoff = 1, qvalueCutoff = 1,
             pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500)
  df <- as.data.frame(res)
  df[order(df$pvalue), ]
}

ontologies <- list(
  DO = function(...) enrichDO(...),
  DGN = function(...) enrichDGN(...),
  KEGG = function(gene, universe, ...) enrichKEGG(gene = gene, universe = universe, organism = "hsa",
                                                  keyType = "ncbi-geneid", ...)
)

summary_rows <- list()
for (ont in names(ontologies)) {
  for (bg_name in names(backgrounds)) {
    df <- run_one(ontologies[[ont]], backgrounds[[bg_name]])
    if ("geneID" %in% names(df)) {
      df$geneSymbol <- vapply(strsplit(df$geneID, "/"), function(ids) {
        paste(gene_table$symbol[match(ids, gene_table$entrez)], collapse = "/")
      }, character(1))
    }
    df$rank <- seq_len(nrow(df))
    out <- file.path(out_dir, sprintf("dose_%s_%s_bg_full.tsv", ont, bg_name))
    write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      ontology = ont, background = bg_name, query_genes_in_ontology = sub(".*/", "", df$GeneRatio[1]),
      background_genes_in_ontology = sub(".*/", "", df$BgRatio[1]), terms_tested = nrow(df),
      terms_fdr_lt_0.05 = sum(df$p.adjust < 0.05), top_term = df$Description[1]
    )
    message(sprintf("%s / %s: %d terms; top = %s", ont, bg_name, nrow(df), df$Description[1]))
  }
}
write.table(do.call(rbind, summary_rows), file.path(out_dir, "dose_run_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("DOSE", packageVersion("DOSE")), paste("clusterProfiler", packageVersion("clusterProfiler")),
             paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db")), R.version.string,
             paste("run_date", Sys.Date())), file.path(out_dir, "dose_versions.txt"))
