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
    top_n = "12",
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

plot_disease_theme_summary <- function(theme_df, title_text, out_prefix) {
  plot_df <- theme_df[theme_df$significant_term_count > 0, ]
  if (nrow(plot_df) == 0) return(invisible(NULL))
  plot_df$best_neg_log10_q <- -log10(pmax(plot_df$best_q_value, .Machine$double.xmin))
  plot_df <- plot_df[order(plot_df$significant_term_fraction), ]
  plot_df$disease_theme <- factor(plot_df$disease_theme, levels = plot_df$disease_theme)
  p <- ggplot(plot_df, aes(x = significant_term_fraction, y = disease_theme)) +
    geom_col(fill = "#1b9e77") +
    geom_point(aes(color = best_neg_log10_q, size = significant_term_count)) +
    scale_color_distiller(palette = "YlOrRd", direction = 1, name = expression(-log[10]("best q"))) +
    geom_text(aes(label = paste0(significant_term_count, "/", matched_term_count, " | ", best_term)), hjust = -0.02, size = 2.8) +
    labs(
      title = title_text,
      x = "Fraction of matched disease terms significant after BH correction",
      y = NULL,
      color = expression(-log[10]("best q")),
      size = "Significant terms"
    ) +
    theme_bw(base_family = "sans") +
    theme(
      plot.title = element_text(size = 13),
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank()
    ) +
    expand_limits(x = min(1.0, max(plot_df$significant_term_fraction, na.rm = TRUE) * 1.2))
  save_plot(p, out_prefix, width = 11, height = 6.8)
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
    theme_bw(base_family = "sans") +
    theme(
      plot.title = element_text(size = 13),
      axis.text.y = element_text(size = 9)
    )
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
          ggtitle(paste0(target_name, ": disease term-gene network"))
        save_plot(fig_cnet, file.path(outdir, paste0(safe_filename(target_name), "_disease_cnetplot")), width = 11, height = 8.5)
      },
      error = function(e) warning(sprintf("Failed cnetplot for %s: %s", target_name, e$message))
    )
    tryCatch(
      {
        fig_heat <- heatplot(enr, showCategory = network_show_n) +
          ggtitle(paste0(target_name, ": disease term-gene heatmap"))
        save_plot(fig_heat, file.path(outdir, paste0(safe_filename(target_name), "_disease_heatplot")), width = 11, height = 8.5)
      },
      error = function(e) warning(sprintf("Failed heatplot for %s: %s", target_name, e$message))
    )
  }

  theme_df <- summarize_disease_themes(result_df, q_cutoff = q_cutoff)
  if (nrow(theme_df) > 0) {
    theme_rows[[target_name]] <- cbind(target = target_name, theme_df)
    write.table(
      theme_df,
      file = file.path(outdir, paste0(safe_filename(target_name), "_disease_theme_summary.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    if (target_name == "All methyl proteins") {
      plot_disease_theme_summary(
        theme_df,
        "All methyl proteins: consolidated disease-theme summary",
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
      geom_col(fill = "#d95f02") +
      geom_text(aes(label = paste0("n=", Count, ", q=", format(p.adjust, digits = 2, scientific = TRUE))), hjust = -0.02, size = 2.8) +
      labs(
        title = paste0(target_name, ": neurodegenerative disease associations"),
        x = "log2 fold-enrichment vs reviewed-human background",
        y = NULL
      ) +
      theme_bw(base_family = "sans") +
      theme(
        plot.title = element_text(size = 13),
        axis.text.y = element_text(size = 9),
        panel.grid.major.y = element_blank()
      ) +
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
