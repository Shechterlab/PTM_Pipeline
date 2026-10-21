#!/usr/bin/env Rscript
# Over-representation tests of protein or gene sets with published tools, for the analyses that ask which sets the
# methylarginine-containing proteins fall in more often than the background (scripts/enrichment.py writes the jobs):
#
#   Rscript scripts/enrich_sets.R JOBS.tsv OUTDIR
#
# JOBS.tsv has one row per job:
#   job            name; the results go to OUTDIR/<job>.tsv
#   tool           enricher: clusterProfiler::enricher on the job's sets (TERM2GENE)
#                  enrichDGN: DOSE::enrichDGN on DOSE's bundled DisGeNET data (the job's sets are not used)
#   target         file of the target identifiers, one per line
#   universe       file of the background identifiers, one per line
#   sets           tab-separated term and gene (TERM2GENE), with a header; for enricher
#   names          tab-separated term and name (TERM2NAME), with a header, or empty
#   whole_universe true: the background is the whole universe. clusterProfiler keeps only the universe identifiers that
#                  some set contains, so the runner adds the universe as one more set (dropped from the results);
#                  false: the background is the universe identifiers in at least one set, clusterProfiler's default
#   min_size, max_size   the set sizes tested (minGSSize, maxGSSize)
#
# Every tested set is written, whatever its p value (pvalueCutoff = qvalueCutoff = 1), with clusterProfiler's columns:
# GeneRatio k/n, BgRatio M/N, FoldEnrichment, RichFactor, the one-sided hypergeometric pvalue, BH p.adjust within the
# job, qvalue and the target identifiers in the set (geneID). OUTDIR/versions.tsv records the R and package versions.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DOSE)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: enrich_sets.R JOBS.tsv OUTDIR")
jobs <- read.delim(args[[1]], colClasses = "character", na.strings = character())
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ids <- function(path) unique(trimws(readLines(path)))
columns <- c("ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment", "RichFactor", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
ALL <- "__whole_universe__"

for (i in seq_len(nrow(jobs))) {
  j <- as.list(jobs[i, ])
  target <- ids(j$target)
  universe <- ids(j$universe)
  min_size <- as.numeric(j$min_size)
  max_size <- as.numeric(j$max_size)
  if (j$tool == "enrichDGN") {
    res <- enrichDGN(target, universe = universe, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                     minGSSize = min_size, maxGSSize = max_size, readable = FALSE)
  } else if (j$tool == "enricher") {
    sets <- read.delim(j$sets, colClasses = "character")[, 1:2]
    names(sets) <- c("term", "gene")
    if (identical(tolower(j$whole_universe), "true")) {
      # the whole-universe set must pass maxGSSize; the job's own sets above max_size are dropped after the test
      sets <- rbind(sets, data.frame(term = ALL, gene = universe))
      max_size <- max(max_size, length(universe))
    }
    term_names <- if (nzchar(j$names)) read.delim(j$names, colClasses = "character")[, 1:2] else NA
    res <- enricher(target, universe = universe, TERM2GENE = sets, TERM2NAME = term_names, pvalueCutoff = 1,
                    qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = min_size, maxGSSize = max_size)
  } else {
    stop(sprintf("job %s: unknown tool %s", j$job, j$tool))
  }
  d <- if (is.null(res)) data.frame() else as.data.frame(res)
  if (nrow(d)) {
    d <- d[d$ID != ALL, , drop = FALSE]
    if (identical(tolower(j$whole_universe), "true")) {
      d <- d[as.numeric(sub("/.*", "", d$BgRatio)) <= as.numeric(j$max_size), , drop = FALSE]
      # BH (stats::p.adjust) over the job's own sets, without the whole-universe set; clusterProfiler's qvalue
      # counted that set too, so it is left out
      d$p.adjust <- p.adjust(d$pvalue, method = "BH")
      d$qvalue <- NA
    }
    d <- d[, intersect(columns, names(d)), drop = FALSE]
  }
  write.table(d, file.path(outdir, paste0(j$job, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("%s: %s, %d sets tested, %d with p.adjust <= 0.05\n", j$job, j$tool, nrow(d), if (nrow(d)) sum(d$p.adjust <= 0.05) else 0L))
}
versions <- data.frame(package = c("R", "clusterProfiler", "DOSE"),
                       version = c(paste(R.version$major, R.version$minor, sep = "."),
                                   as.character(packageVersion("clusterProfiler")), as.character(packageVersion("DOSE"))))
write.table(versions, file.path(outdir, "versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
