#!/usr/bin/env Rscript
# Export, once, the Disease Ontology data DOSE's enrichDO uses, so that the build tests the same gene sets without
# downloading them (DOSE 4.4 fetches them at run time, through GOSemSim, from YuLab's HDO.sqlite, which follows the
# Human Disease Ontology releases):
#
#   hdo_<date>_term_genes.tsv   one row per Disease Ontology term and Entrez gene: the genes of the term and of every
#                               term under it (DOSE's term-to-all-genes table, what enrichDO tests)
#   hdo_<date>_terms.tsv        every term with its name and its parents (is_a)
#   hdo_<date>_source.tsv       the ontology release, its source and the downloaded file's MD5
#
#   Rscript scripts/00_export_disease_ontology.R inputs/disease
#
# <date> is the Disease Ontology release in the file (SOURCEDATE), e.g. 2026-09-30.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: 00_export_disease_ontology.R OUTDIR")
outdir <- args[[1]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(DOSE)
  library(GOSemSim)
  library(DBI)
  library(AnnotationDbi)
})

term_genes <- DOSE:::get_ont2allgene("HDO")
db <- GOSemSim:::load_onto("HDO")
con <- dbconn(db)
meta <- dbReadTable(con, "metadata")
release <- meta$value[meta$name == "SOURCEDATE"]
date <- sprintf("%s-%s-%s", substr(release, 1, 4), substr(release, 5, 6), substr(release, 7, 8))
terms <- dbReadTable(con, "term")
parents <- dbReadTable(con, "parent")
names(terms)[1:2] <- c("term_id", "name")
names(parents)[1:2] <- c("term_id", "parent_id")
parent_of <- tapply(parents$parent_id, parents$term_id, function(x) paste(sort(unique(x)), collapse = ";"))
terms$parents <- unname(ifelse(is.na(parent_of[terms$term_id]), "", parent_of[terms$term_id]))
terms <- terms[order(terms$term_id), c("term_id", "name", "parents")]

genes <- data.frame(term_id = rep(names(term_genes), lengths(term_genes)), entrez_id = unlist(term_genes, use.names = FALSE),
                    stringsAsFactors = FALSE)
genes <- unique(genes[!is.na(genes$entrez_id) & genes$entrez_id != "", ])
genes <- genes[order(genes$term_id, as.integer(genes$entrez_id)), ]

file <- file.path(yulab.utils::user_dir("GOSemSim"), "HDO.sqlite")
source <- data.frame(
  field = c("ontology", "release", "source_url", "database", "database_md5", "DOSE", "GOSemSim"),
  value = c(meta$value[meta$name == "SOURCENAME"], date, meta$value[meta$name == "SOURCURL"],
            "HDO.sqlite (YuLab, https://yulab-smu.top/DOSE), downloaded by GOSemSim", unname(tools::md5sum(file)),
            as.character(packageVersion("DOSE")), as.character(packageVersion("GOSemSim"))))
write.table(genes, file.path(outdir, sprintf("hdo_%s_term_genes.tsv", date)), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(terms, file.path(outdir, sprintf("hdo_%s_terms.tsv", date)), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(source, file.path(outdir, sprintf("hdo_%s_source.tsv", date)), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Disease Ontology", date, ":", length(unique(genes$term_id)), "terms with genes,", nrow(genes), "term-gene rows,",
    nrow(terms), "terms\n")
