#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(ggplot2)
})

parse_args <- function(args) {
  defaults <- list(
    labels = "PTM_results/functional_class_union/protein_functional_labels.tsv",
    outdir = "PTM_results/functional_class_union/clusterprofiler_go",
    label_column = "functional_class_broad",
    targets = "Other / unclassified,Poorly characterized / specialized",
    include_all_methyl = "false",
    ontologies = "BP,CC,MF",
    min_size = "10",
    max_size = "5000",
    top_n = "6",
    simplify_cutoff = "0.35",
    padj_cutoff = "0.05"
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

extract_primary_symbol <- function(symbol_col, genename_col) {
  values <- trimws(ifelse(is.na(symbol_col), "", symbol_col))
  extracted <- rep("", length(values))

  nonempty_symbol <- nzchar(values)
  extracted[nonempty_symbol] <- sub("\\s*;.*$", "", values[nonempty_symbol])

  need_fallback <- !nzchar(extracted)
  if (any(need_fallback)) {
    fallback <- trimws(ifelse(is.na(genename_col), "", genename_col))
    fallback <- sub("^.*?Name:\\s*", "", fallback)
    fallback <- sub(";.*$", "", fallback)
    extracted[need_fallback] <- fallback[need_fallback]
  }

  extracted <- trimws(extracted)
  extracted[grepl("\\s", extracted)] <- ""
  extracted
}

split_labels <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

has_label <- function(label_string, target_label) {
  target_label %in% split_labels(label_string)
}

safe_filename <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", x)
  y <- gsub("^_+|_+$", "", y)
  tolower(y)
}

save_plot <- function(plot_obj, out_prefix, width = 10, height = 7) {
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

make_plot <- function(df, title_text, out_prefix, top_n) {
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }

  top_df <- do.call(
    rbind,
    lapply(split(df, df$ONTOLOGY), function(chunk) {
      chunk <- chunk[order(chunk$p.adjust, -chunk$Count), ]
      head(chunk, top_n)
    })
  )
  rownames(top_df) <- NULL

  top_df$neg_log10_padj <- -log10(pmax(top_df$p.adjust, .Machine$double.xmin))
  top_df$plot_label <- paste(top_df$Description, sprintf("(%s)", top_df$ONTOLOGY))
  top_df <- top_df[order(top_df$ONTOLOGY, top_df$neg_log10_padj), ]
  top_df$plot_label <- factor(top_df$plot_label, levels = top_df$plot_label)

  p <- ggplot(top_df, aes(x = neg_log10_padj, y = plot_label, size = Count, color = ONTOLOGY)) +
    geom_point(alpha = 0.9) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      title = title_text,
      x = expression(-log[10]("BH-adjusted p")),
      y = NULL,
      size = "Gene count"
    ) +
    theme_minimal(base_family = "sans", base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  save_plot(p, out_prefix, width = 10, height = 7)
}

make_emap_plot <- function(ego, title_text, out_prefix, show_n) {
  result_df <- as.data.frame(ego)
  if (nrow(result_df) < 4) {
    return(invisible(NULL))
  }

  result_df <- result_df[order(result_df$p.adjust, -result_df$Count), ]
  result_df <- head(result_df, show_n)
  if (nrow(result_df) < 4) {
    return(invisible(NULL))
  }

  ego_top <- ego
  ego_top@result <- ego_top@result[match(result_df$ID, ego_top@result$ID), , drop = FALSE]

  ego_simple <- tryCatch(
    simplify(ego_top, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
  sim_obj <- tryCatch(
    pairwise_termsim(ego_simple, showCategory = min(show_n, nrow(as.data.frame(ego_simple)))),
    error = function(e) NULL
  )
  if (is.null(sim_obj) || nrow(as.data.frame(sim_obj)) < 4) {
    return(invisible(NULL))
  }
  p <- emapplot(
    sim_obj,
    showCategory = min(show_n, nrow(as.data.frame(sim_obj))),
    node_label = "group",
    layout = "kk"
  ) +
    ggtitle(title_text)
  save_plot(p, out_prefix, width = 11, height = 8.5)
}

simplify_for_plot <- function(ego, ont, cutoff) {
  result_df <- as.data.frame(ego)
  if (nrow(result_df) == 0) {
    return(result_df)
  }
  simplified <- tryCatch(
    simplify(ego, cutoff = cutoff, by = "p.adjust", select_fun = min),
    error = function(e) ego
  )
  out <- as.data.frame(simplified)
  if (nrow(out) == 0) {
    out <- result_df
  }
  out <- out[out$p.adjust <= padj_cutoff, , drop = FALSE]
  out$ONTOLOGY <- ont
  out[order(out$p.adjust, -out$Count), ]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
labels_path <- args$labels
outdir <- args$outdir
label_column <- args$label_column
targets <- trimws(unlist(strsplit(args$targets, ",", fixed = TRUE)))
include_all_methyl <- tolower(args$include_all_methyl) %in% c("1", "true", "yes", "y")
ontologies <- trimws(unlist(strsplit(args$ontologies, ",", fixed = TRUE)))
min_size <- as.integer(args$min_size)
max_size <- as.integer(args$max_size)
top_n <- as.integer(args$top_n)
simplify_cutoff <- as.numeric(args$simplify_cutoff)
padj_cutoff <- as.numeric(args$padj_cutoff)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

labels_df <- read.delim(labels_path, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

required_cols <- c("canonical_UniProtAC", "substrate_genename", "genename", "is_arg_methyl_protein", label_column)
missing_cols <- setdiff(required_cols, colnames(labels_df))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

labels_df$primary_symbol <- extract_primary_symbol(labels_df$substrate_genename, labels_df$genename)
labels_df <- labels_df[nzchar(labels_df$primary_symbol), ]
labels_df$is_arg_methyl_protein <- as.logical(labels_df$is_arg_methyl_protein)

unique_symbols <- unique(labels_df$primary_symbol)
symbol_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique_symbols,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
  )
)
symbol_map <- symbol_map[!is.na(symbol_map$ENTREZID) & !is.na(symbol_map$SYMBOL), c("SYMBOL", "ENTREZID")]
symbol_map <- unique(symbol_map)

labels_df <- merge(labels_df, symbol_map, by.x = "primary_symbol", by.y = "SYMBOL", all.x = TRUE)
universe_entrez <- unique(labels_df$ENTREZID[!is.na(labels_df$ENTREZID)])

summary_rows <- list()
target_specs <- lapply(targets, function(x) list(name = x, mode = "label"))
if (include_all_methyl) {
  target_specs <- c(list(list(name = "All methyl proteins", mode = "all_methyl")), target_specs)
}

for (spec in target_specs) {
  target <- spec$name
  message(sprintf("Processing target: %s", target))
  methyl_mask <- labels_df$is_arg_methyl_protein
  if (identical(spec$mode, "all_methyl")) {
    target_df <- labels_df[methyl_mask, ]
  } else {
    target_mask <- vapply(labels_df[[label_column]], has_label, logical(1), target_label = target)
    target_df <- labels_df[target_mask & methyl_mask, ]
  }

  target_symbols <- unique(target_df$primary_symbol)
  target_entrez <- unique(target_df$ENTREZID[!is.na(target_df$ENTREZID)])

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    target_label = target,
    universe_proteins = length(unique(labels_df$canonical_UniProtAC)),
    universe_symbols = length(unique(labels_df$primary_symbol)),
    universe_entrez = length(universe_entrez),
    target_proteins = length(unique(target_df$canonical_UniProtAC)),
    target_symbols = length(target_symbols),
    target_entrez = length(target_entrez),
    stringsAsFactors = FALSE
  )

  safe_target <- safe_filename(target)
  message(sprintf("  proteins=%d symbols=%d entrez=%d", length(unique(target_df$canonical_UniProtAC)), length(target_symbols), length(target_entrez)))
  write.table(
    unique(target_df[, c("canonical_UniProtAC", "primary_symbol", "protein_name", label_column, "ENTREZID")]),
    file = file.path(outdir, paste0(safe_target, "_input_genes.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  if (length(target_entrez) < min_size) {
    write.table(
      data.frame(),
      file = file.path(outdir, paste0(safe_target, "_go_combined.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    next
  }

  combined_results <- list()
  combined_plot_results <- list()
  for (ont in ontologies) {
    message(sprintf("  enrichGO ont=%s", ont))
    ego <- suppressMessages(
      enrichGO(
        gene = target_entrez,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = min_size,
        maxGSSize = max_size,
        readable = TRUE
      )
    )

    result_df <- as.data.frame(ego)
    message(sprintf("    rows=%d", nrow(result_df)))
    if (nrow(result_df) == 0) {
      write.table(
        data.frame(),
        file = file.path(outdir, paste0(safe_target, "_", tolower(ont), ".tsv")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
      next
    }

    result_df <- result_df[order(result_df$p.adjust, -result_df$Count), ]
    write.table(
      result_df,
      file = file.path(outdir, paste0(safe_target, "_", tolower(ont), ".tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    result_df$ONTOLOGY <- ont
    combined_results[[length(combined_results) + 1]] <- result_df
    combined_plot_results[[length(combined_plot_results) + 1]] <- simplify_for_plot(ego, ont, simplify_cutoff)

    if (ont == "BP") {
      make_emap_plot(
        ego,
        sprintf("GO BP term network: %s", target),
        file.path(outdir, paste0(safe_target, "_bp_emapplot")),
        show_n = max(top_n, 10)
      )
    }
  }

  if (length(combined_results) == 0) {
    write.table(
      data.frame(),
      file = file.path(outdir, paste0(safe_target, "_go_combined.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    next
  }

  combined_df <- do.call(rbind, combined_results)
  rownames(combined_df) <- NULL
  write.table(
    combined_df,
    file = file.path(outdir, paste0(safe_target, "_go_combined.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  combined_plot_results <- combined_plot_results[vapply(combined_plot_results, nrow, integer(1)) > 0]
  combined_plot_df <- if (length(combined_plot_results) > 0) do.call(rbind, combined_plot_results) else data.frame()
  rownames(combined_plot_df) <- NULL
  write.table(
    combined_plot_df,
    file = file.path(outdir, paste0(safe_target, "_go_combined_simplified.tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  if (nrow(combined_plot_df) > 0) {
    make_plot(
      combined_plot_df,
      sprintf("clusterProfiler GO Enrichment: %s", target),
      file.path(outdir, paste0(safe_target, "_go_dotplot")),
      top_n = top_n
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
write.table(
  summary_df,
  file = file.path(outdir, "clusterprofiler_input_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
