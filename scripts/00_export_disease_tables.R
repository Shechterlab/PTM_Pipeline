#!/usr/bin/env Rscript
# Export, once, the tables the disease enrichment needs, so the pipeline itself runs without R:
#
#   dose_<version>_disgenet_term_genes.tsv   DOSE's DisGeNET data: one row per term and Entrez gene
#                                            (DGN_PATHID2EXTID, with DGN_PATHID2NAME for the names)
#   orghsegdb_<version>_symbol_entrez.tsv    org.Hs.eg.db: every gene symbol with its Entrez gene, and the symbol
#                                            org.Hs.eg.db gives each Entrez gene (used to name genes)
#
#   Rscript scripts/00_export_disease_tables.R inputs/disease
#
# DOSE's enrichDGN reads genes per term from DGN_PATHID2EXTID and terms per gene from DGN_EXTID2PATHID; the export
# stops unless the two are exact inverses, so one table serves both.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 00_export_disease_tables.R OUTDIR")
outdir <- args[[1]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(DOSE)
  library(org.Hs.eg.db)
})

data("DGN_PATHID2EXTID", package = "DOSE", envir = environment())
data("DGN_EXTID2PATHID", package = "DOSE", envir = environment())
data("DGN_PATHID2NAME", package = "DOSE", envir = environment())
term_ids <- names(DGN_PATHID2EXTID)
terms <- data.frame(
  term_id = rep(term_ids, lengths(DGN_PATHID2EXTID)),
  term_name = unname(DGN_PATHID2NAME[rep(term_ids, lengths(DGN_PATHID2EXTID))]),
  entrez_id = unlist(DGN_PATHID2EXTID, use.names = FALSE),
  stringsAsFactors = FALSE
)
terms <- unique(terms)
inverse <- unique(data.frame(
  entrez_id = rep(names(DGN_EXTID2PATHID), lengths(DGN_EXTID2PATHID)),
  term_id = unlist(DGN_EXTID2PATHID, use.names = FALSE),
  stringsAsFactors = FALSE
))
a <- paste(terms$term_id, terms$entrez_id)
b <- paste(inverse$term_id, inverse$entrez_id)
if (!setequal(a, b)) stop("DGN_PATHID2EXTID and DGN_EXTID2PATHID are not inverses")
terms <- terms[order(terms$term_id, as.integer(terms$entrez_id)), ]
dose <- as.character(packageVersion("DOSE"))
write.table(terms, file.path(outdir, paste0("dose_", dose, "_disgenet_term_genes.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
                                 keytype = "SYMBOL", columns = "ENTREZID")
symbols <- unique(symbols[!is.na(symbols$ENTREZID), c("SYMBOL", "ENTREZID")])
names(symbols) <- c("symbol", "entrez_id")
# The symbol org.Hs.eg.db gives an Entrez gene (what setReadable shows), marked per row.
named <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(terms$entrez_id), keytype = "ENTREZID", columns = "SYMBOL")
named <- named[!is.na(named$SYMBOL) & !duplicated(named$ENTREZID), ]
symbols$entrez_symbol <- ifelse(paste(symbols$entrez_id, symbols$symbol) %in% paste(named$ENTREZID, named$SYMBOL), "true", "false")
symbols <- symbols[order(symbols$symbol, as.integer(symbols$entrez_id)), ]
org <- as.character(packageVersion("org.Hs.eg.db"))
write.table(symbols, file.path(outdir, paste0("orghsegdb_", org, "_symbol_entrez.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DOSE", dose, ":", length(unique(terms$term_id)), "terms,", nrow(terms), "term-gene rows\n")
cat("org.Hs.eg.db", org, ":", nrow(symbols), "symbol-gene rows\n")
