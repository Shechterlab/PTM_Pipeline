#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 34_export_dgn_mapping.R OUTPUT_TSV")
suppressPackageStartupMessages(library(DOSE))
data("DGN_PATHID2EXTID", package = "DOSE", envir = environment())
data("DGN_PATHID2NAME", package = "DOSE", envir = environment())
term_ids <- intersect(names(DGN_PATHID2NAME), names(DGN_PATHID2EXTID))
out <- data.frame(
  term_id = rep(term_ids, lengths(DGN_PATHID2EXTID[term_ids])),
  term_name = rep(unname(DGN_PATHID2NAME[term_ids]), lengths(DGN_PATHID2EXTID[term_ids])),
  entrez_id = unlist(DGN_PATHID2EXTID[term_ids], use.names = FALSE),
  stringsAsFactors = FALSE
)
write.table(out, args[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
