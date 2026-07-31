#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DOSE)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
})

parse_args <- function(args) {
  defaults <- list(
    labels = "",
    outdir = "PTM_results/functional_class_union/disease_enrichment",
    q_cutoff = "0.05",
    top_n = "8",
    min_size = "10",
    max_size = "300"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    name <- substring(key, 3)
    if (!(name %in% names(defaults))) {
      stop(sprintf("Unknown argument: %s", key))
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for %s", key))
    }
    defaults[[name]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

safe_filename <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  tolower(y)
}

save_plot <- function(plot_obj, out_prefix, width = 9, height = 7) {
  ggsave(paste0(out_prefix, ".png"), plot_obj, width = width, height = height, dpi = 300)
  if (tolower(Sys.getenv("PTM_SKIP_PDF", "0")) %in% c("1", "true", "yes", "y")) {
    return(invisible(NULL))
  }
  pdf_device <- if (capabilities("cairo")) cairo_pdf else grDevices::pdf
  tryCatch(
    ggsave(paste0(out_prefix, ".pdf"), plot_obj, width = width, height = height, device = pdf_device),
    error = function(e) warning(sprintf("Failed to write PDF %s: %s", out_prefix, e$message))
  )
}

theme_ptm <- function() {
  theme_classic(base_family = "Arial", base_size = 11) +
    theme(
      plot.title = element_text(size = 12.5, family = "Arial"),
      axis.title = element_text(size = 11, family = "Arial"),
      axis.text = element_text(size = 10, family = "Arial", colour = "#111111"),
      legend.title = element_text(size = 9.5, family = "Arial"),
      legend.text = element_text(size = 9, family = "Arial"),
      panel.grid = element_blank()
    )
}

format_compact_q <- function(x) {
  if (is.na(x)) return("q=NA")
  if (x == 0) return("q≈0")
  if (x < 1e-2) {
    exponent <- floor(log10(x))
    mantissa <- round(x / (10 ^ exponent), 0)
    if (mantissa >= 10) {
      mantissa <- mantissa / 10
      exponent <- exponent + 1
    }
    superscript_map <- c(
      "-" = "\u207b",
      "0" = "\u2070",
      "1" = "\u00b9",
      "2" = "\u00b2",
      "3" = "\u00b3",
      "4" = "\u2074",
      "5" = "\u2075",
      "6" = "\u2076",
      "7" = "\u2077",
      "8" = "\u2078",
      "9" = "\u2079"
    )
    exponent_chars <- strsplit(as.character(exponent), "", fixed = TRUE)[[1]]
    exponent_label <- paste0(unname(superscript_map[exponent_chars]), collapse = "")
    return(paste0("q=", mantissa, "\u00D710", exponent_label))
  }
  if (x < 0.1) return(sprintf("q=%.2f", x))
  if (x < 1) return(sprintf("q=%.1f", x))
  sprintf("q=%.0f", x)
}

format_q_threshold_note <- function(x) {
  if (!is.na(x) && x <= 0.05) return("All shown BH q<0.05")
  label <- format_compact_q(x)
  sub("^q=", "All shown q\u2264", label)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$labels) || args$labels == "") {
  stop("--labels is required")
}
q_cutoff <- as.numeric(args$q_cutoff)
top_n <- as.integer(args$top_n)
min_size <- as.integer(args$min_size)
max_size <- as.integer(args$max_size)

primary_gene <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  strsplit(as.character(x), ";", fixed = TRUE)[[1]][1]
}

parse_ratio_num <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "/", fixed = TRUE)
  sapply(parts, function(v) {
    if (length(v) != 2) return(NA_real_)
    as.numeric(v[1]) / as.numeric(v[2])
  })
}

add_fold_enrichment <- function(df) {
  if (nrow(df) == 0) return(df)
  df$gene_ratio_num <- parse_ratio_num(df$GeneRatio)
  df$bg_ratio_num <- parse_ratio_num(df$BgRatio)
  df$fold_enrichment <- df$gene_ratio_num / df$bg_ratio_num
  df$log2_fold_enrichment <- log2(df$fold_enrichment)
  df
}

to_bool <- function(x) {
  values <- trimws(tolower(as.character(x)))
  values %in% c("1", "true", "t", "yes", "y")
}

neuro_pattern <- paste(
  c(
    "amyotrophic", "als", "motor neuron", "frontotemporal", "ftd",
    "neurodeg", "alzheimer", "parkinson", "dementia", "ataxia",
    "huntington", "tauopath", "inclusion body", "spinocerebellar",
    "neuron degeneration", "neurofilament"
  ),
  collapse = "|"
)

disease_theme_patterns <- list(
  "ALS / motor neuron / FTD" = "amyotrophic|als|motor neuron|frontotemporal|ftd|bulbar palsy",
  "Alzheimer / tau / dementia" = "alzheimer|tauopath|dementia|semantic dementia|cognitive disorder",
  "Ataxia / cerebellar degeneration" = "ataxia|spinocerebellar|friedreich|cerebellar",
  "Parkinsonism" = "parkinson",
  "Cancer / neoplasia" = "cancer|carcinoma|sarcoma|glioma|leukemia|lymphoma|melanoma|tumou?r|adenoma|neoplas|myeloma",
  "Myopathy / muscle / cardiomyopathy" = "myopath|myopathy|muscular dystrophy|cardiomyopathy|rhabdomyolysis|myositis",
  "Developmental / congenital syndromes" = "congenital|dysmorph|cranio|palate|facies|syndrom|hypoplasia|epicanthus|forehead|midface|feeding difficulties",
  "Premature aging / genome instability" = "werner|cockayne|ataxia telangiectasia|xeroderma|fanconi|bloom|genome instability",
  "Immune / inflammatory" = "autoimmune|inflamm|lupus|arthritis|colitis|crohn|asthma|dermatitis|psoriasis",
  "Metabolic / mitochondrial" = "mitochond|metabolic|glycogen|lipodystrophy|diabetes|obesity"
)

summarize_disease_themes <- function(df, q_cutoff) {
  if (nrow(df) == 0) return(data.frame())
  sig_df <- df[df$p.adjust <= q_cutoff, ]
  rows <- list()
  for (theme_name in names(disease_theme_patterns)) {
    pattern <- disease_theme_patterns[[theme_name]]
    all_matches <- df[grepl(pattern, df$Description, ignore.case = TRUE), ]
    sig_matches <- sig_df[grepl(pattern, sig_df$Description, ignore.case = TRUE), ]
    if (nrow(all_matches) == 0) next
    top_match <- all_matches[order(all_matches$p.adjust, -all_matches$Count), ][1, , drop = FALSE]
    rows[[length(rows) + 1]] <- data.frame(
      disease_theme = theme_name,
      matched_term_count = nrow(all_matches),
      significant_term_count = nrow(sig_matches),
      significant_term_fraction = nrow(sig_matches) / nrow(all_matches),
      best_term = top_match$Description[[1]],
      best_q_value = top_match$p.adjust[[1]],
      best_log2_fold_enrichment = top_match$log2_fold_enrichment[[1]],
      median_log2_fold_enrichment = if (nrow(sig_matches) > 0) median(sig_matches$log2_fold_enrichment, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  out <- out[order(-out$significant_term_count, out$best_q_value), ]
  rownames(out) <- NULL
  out
}

fisher_odds <- function(a, b, c, d) {
  ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
}

theme_union_enrichment <- function(target_entrez, universe_entrez) {
  get_dgn <- getFromNamespace("get_DGN_data", "DOSE")
  dgn_env <- get_dgn()
  pathid2name <- get("PATHID2NAME", envir = dgn_env)
  pathid2extid <- get("PATHID2EXTID", envir = dgn_env)

  target_set <- unique(as.character(target_entrez))
  universe_set <- unique(as.character(universe_entrez))
  background_set <- setdiff(universe_set, target_set)

  rows <- list()
  for (theme_name in names(disease_theme_patterns)) {
    pattern <- disease_theme_patterns[[theme_name]]
    matched_ids <- names(pathid2name)[grepl(pattern, unname(pathid2name), ignore.case = TRUE)]
    if (length(matched_ids) == 0) next
    theme_genes <- unique(unlist(pathid2extid[matched_ids], use.names = FALSE))
    theme_genes <- intersect(as.character(theme_genes), universe_set)
    if (length(theme_genes) == 0) next
    a <- sum(target_set %in% theme_genes)
    b <- length(target_set) - a
    c <- sum(background_set %in% theme_genes)
    d <- length(background_set) - c
    rows[[length(rows) + 1]] <- data.frame(
      disease_theme = theme_name,
      matched_dgn_term_count = length(matched_ids),
      target_count = a,
      background_count = c,
      total_theme_gene_count = a + c,
      target_fraction_in_theme = if ((a + c) > 0) a / (a + c) else NA_real_,
      odds_ratio = fisher_odds(a, b, c, d),
      log2_odds_ratio = log2(fisher_odds(a, b, c, d)),
      p_value = fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "two.sided")$p.value,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out <- out[out$target_count > 0, ]
  out <- out[order(out$odds_ratio, out$target_count), ]
  rownames(out) <- NULL
  out
}

plot_disease_theme_summary <- function(theme_df, title_text, out_prefix) {
  plot_df <- theme_df[theme_df$q_value <= q_cutoff, ]
  if (nrow(plot_df) == 0) return(invisible(NULL))
  plot_df <- plot_df[order(plot_df$odds_ratio), ]
  plot_df$disease_theme <- factor(plot_df$disease_theme, levels = plot_df$disease_theme)
  max_x <- max(plot_df$odds_ratio, na.rm = TRUE)
  max_q <- max(plot_df$q_value, na.rm = TRUE)
  p <- ggplot(plot_df, aes(y = disease_theme)) +
    geom_col(aes(x = odds_ratio), fill = "#8c8c8c", width = 0.72) +
    geom_text(
      aes(x = odds_ratio + max_x * 0.03, label = paste0(target_count, "/", total_theme_gene_count)),
      hjust = 0,
      size = 3.3,
      family = "Arial"
    ) +
    labs(
      title = title_text,
      x = "Odds ratio vs reviewed-human background",
      y = NULL
    ) +
    theme_ptm() +
    coord_cartesian(xlim = c(0, max_x * 1.22), clip = "off") +
    annotate("text", x = max_x * 1.21, y = 0.55, label = format_q_threshold_note(max_q), hjust = 1, vjust = 0, size = 3.0, family = "Arial")
  save_plot(p, out_prefix, width = 10.8, height = 6.4)
}

labels_df <- read.delim(args$labels, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
labels_df$is_arg_methyl_protein <- to_bool(labels_df$is_arg_methyl_protein)
labels_df$gene_symbol <- vapply(labels_df$substrate_genename, primary_gene, FUN.VALUE = character(1))
missing_symbol <- is.na(labels_df$gene_symbol) | labels_df$gene_symbol == ""
labels_df$gene_symbol[missing_symbol] <- vapply(labels_df$genename[missing_symbol], primary_gene, FUN.VALUE = character(1))
labels_df$gene_symbol <- trimws(labels_df$gene_symbol)
labels_df <- labels_df[!is.na(labels_df$gene_symbol) & labels_df$gene_symbol != "", ]

universe_symbols <- unique(labels_df$gene_symbol)
symbol_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = universe_symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)
symbol_map <- unique(symbol_map[!is.na(symbol_map$ENTREZID), c("SYMBOL", "ENTREZID")])
universe_entrez <- unique(symbol_map$ENTREZID)

targets <- list(
  "All methyl proteins" = labels_df[labels_df$is_arg_methyl_protein, ],
  "Other / unclassified" = labels_df[labels_df$is_arg_methyl_protein & labels_df$functional_class_broad == "Other / unclassified", ],
  "Poorly characterized / specialized" = labels_df[labels_df$is_arg_methyl_protein & labels_df$functional_class_broad == "Poorly characterized / specialized", ]
)

outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

input_summary <- data.frame(
  target = c("Universe", names(targets)),
  protein_count = c(length(unique(labels_df$canonical_UniProtAC)), sapply(targets, function(x) length(unique(x$canonical_UniProtAC)))),
  symbol_count = c(length(universe_symbols), sapply(targets, function(x) length(unique(x$gene_symbol)))),
  entrez_count = c(length(universe_entrez), rep(NA_integer_, length(targets))),
  stringsAsFactors = FALSE
)

combined_rows <- list()
neuro_rows <- list()
theme_rows <- list()

for (target_name in names(targets)) {
  message("Processing disease target: ", target_name)
  target_df <- targets[[target_name]]
  symbols <- unique(target_df$gene_symbol)
  entrez <- unique(symbol_map$ENTREZID[symbol_map$SYMBOL %in% symbols])
  input_summary$entrez_count[input_summary$target == target_name] <- length(entrez)
  if (length(entrez) < 10) {
    next
  }

  enr <- suppressMessages(
    enrichDGN(
      gene = entrez,
      universe = universe_entrez,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = q_cutoff,
      qvalueCutoff = q_cutoff,
      pAdjustMethod = "BH",
      readable = TRUE
    )
  )
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
    next
  }
  result_df <- add_fold_enrichment(as.data.frame(enr))
  result_df$target <- target_name
  result_df <- result_df[order(result_df$p.adjust, -result_df$Count), ]
  write.table(
    result_df,
    file = file.path(outdir, paste0(safe_filename(target_name), "_disease_enrichment.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  combined_rows[[target_name]] <- result_df

  show_n <- min(top_n, nrow(result_df))
  fig <- dotplot(enr, showCategory = show_n) +
    ggtitle(paste0(target_name, ": disease enrichment")) +
    theme_ptm()
  save_plot(fig, file.path(outdir, paste0(safe_filename(target_name), "_disease_dotplot")), width = 9, height = 7)

  filtered_df <- result_df[result_df$p.adjust <= q_cutoff, ]
  if (nrow(filtered_df) >= 3) {
    network_show_n <- min(8, nrow(filtered_df))
    tryCatch(
      {
        fig_cnet <- tryCatch(
          cnetplot(enr, showCategory = network_show_n, colorEdge = TRUE, circular = FALSE, node_label = "all"),
          error = function(e) cnetplot(enr, showCategory = network_show_n, node_label = "all")
        ) +
          ggtitle(paste0(target_name, ": disease term-gene network")) +
          theme_ptm()
        save_plot(fig_cnet, file.path(outdir, paste0(safe_filename(target_name), "_disease_cnetplot")), width = 11, height = 8.5)
      },
      error = function(e) warning(sprintf("Failed cnetplot for %s: %s", target_name, e$message))
    )
    tryCatch(
      {
        fig_heat <- heatplot(enr, showCategory = network_show_n) +
          ggtitle(paste0(target_name, ": disease term-gene heatmap")) +
          theme_ptm()
        save_plot(fig_heat, file.path(outdir, paste0(safe_filename(target_name), "_disease_heatplot")), width = 11, height = 8.5)
      },
      error = function(e) warning(sprintf("Failed heatplot for %s: %s", target_name, e$message))
    )
  }

  theme_term_df <- summarize_disease_themes(result_df, q_cutoff = q_cutoff)
  theme_df <- theme_union_enrichment(entrez, universe_entrez)
  if (nrow(theme_df) > 0) {
    theme_rows[[target_name]] <- cbind(target = target_name, theme_df)
    write.table(
      theme_df,
      file = file.path(outdir, paste0(safe_filename(target_name), "_disease_theme_summary.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    if (nrow(theme_term_df) > 0) {
      write.table(
        theme_term_df,
        file = file.path(outdir, paste0(safe_filename(target_name), "_disease_theme_term_summary.tsv")),
        sep = "\t", row.names = FALSE, quote = FALSE
      )
    }
    if (target_name == "All methyl proteins") {
      plot_disease_theme_summary(
        theme_df,
        "All methyl proteins: consolidated disease-theme enrichment",
        file.path(outdir, paste0(safe_filename(target_name), "_disease_theme_summary"))
      )
    }
  }

  neuro_df <- result_df[grepl(neuro_pattern, result_df$Description, ignore.case = TRUE), ]
  if (nrow(neuro_df) > 0) {
    neuro_df <- neuro_df[seq_len(min(12, nrow(neuro_df))), ]
    neuro_rows[[target_name]] <- neuro_df
    neuro_df$Description <- factor(neuro_df$Description, levels = rev(neuro_df$Description))
    p <- ggplot(neuro_df, aes(x = log2_fold_enrichment, y = Description)) +
      geom_col(fill = "#8c8c8c") +
      geom_text(aes(label = vapply(p.adjust, format_compact_q, character(1))), hjust = -0.02, size = 3.0) +
      labs(
        title = paste0(target_name, ": neurodegenerative disease associations"),
        x = "log2 fold-enrichment vs reviewed-human background",
        y = NULL
      ) +
      theme_ptm() +
      expand_limits(x = max(neuro_df$log2_fold_enrichment, na.rm = TRUE) * 1.25)
    save_plot(p, file.path(outdir, paste0(safe_filename(target_name), "_neurodegenerative_focus")), width = 9, height = 5.8)
    write.table(
      neuro_df,
      file = file.path(outdir, paste0(safe_filename(target_name), "_neurodegenerative_focus.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}

write.table(input_summary, file = file.path(outdir, "disease_input_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (length(combined_rows) > 0) {
  write.table(do.call(rbind, combined_rows), file = file.path(outdir, "disease_enrichment_combined.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
if (length(neuro_rows) > 0) {
  write.table(do.call(rbind, neuro_rows), file = file.path(outdir, "neurodegenerative_disease_focus_combined.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
if (length(theme_rows) > 0) {
  write.table(do.call(rbind, theme_rows), file = file.path(outdir, "disease_theme_summary_combined.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
