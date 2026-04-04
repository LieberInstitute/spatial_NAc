library(spatialLIBD)
library(tidyverse)
library(ComplexHeatmap)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(patchwork)

dat_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration.rds"))

sce_dir <- here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(sce_dir, "sce_CellType_noresiduals.Rds"))

pairwise <- sn_registration$pairwise

# =========================================================
# Subset to DRD1 neurons
# =========================================================
plot_celltypes <- c("DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D")

sce <- sce[, sce$CellType.Final %in% plot_celltypes]
sce$CellType.Final <- factor(as.character(sce$CellType.Final), levels = plot_celltypes)

rownames(sce) <- rowData(sce)$gene_name

genes_to_plot <- c("DRD1", "PDYN", "TSHZ1")

# Use existing cluster colors
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)
celltype_colors <- cluster_cols[levels(sce$CellType.Final)]

# =========================================================
# Settings
# =========================================================
sig_prefix <- "fdr_"
alpha_thresh <- 0.05

outdir <- here("plots", "25_revisions", "DRD1_violin_plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Helper functions
# =========================================================
p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ NA_character_
  )
}

get_expr_mat <- function(sce) {
  if ("logcounts" %in% SummarizedExperiment::assayNames(sce)) {
    SummarizedExperiment::assay(sce, "logcounts")
  } else if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
    message("Using counts assay because logcounts was not found.")
    SummarizedExperiment::assay(sce, "counts")
  } else {
    stop("Neither 'logcounts' nor 'counts' assay found in sce.")
  }
}

get_gene_column <- function(df) {
  if ("gene" %in% colnames(df)) return("gene")
  if ("genes" %in% colnames(df)) return("genes")
  NULL
}

get_gene_expr_df <- function(sce, gene_symbol) {
  expr <- get_expr_mat(sce)

  if (!gene_symbol %in% rownames(expr)) {
    stop(sprintf("Gene '%s' not found in rownames(expr).", gene_symbol))
  }

  tibble(
    cell_type = sce$CellType.Final,
    expr = as.numeric(expr[gene_symbol, ])
  ) %>%
    dplyr::filter(cell_type %in% plot_celltypes) %>%
    dplyr::mutate(
      cell_type = factor(as.character(cell_type), levels = plot_celltypes)
    )
}

get_sig_pairs_for_gene <- function(pairwise, gene_symbol) {
  gene_col <- get_gene_column(pairwise)

  if (!is.null(gene_col)) {
    gene_idx <- which(pairwise[[gene_col]] == gene_symbol)[1]
  } else if (!is.null(rownames(pairwise)) && gene_symbol %in% rownames(pairwise)) {
    gene_idx <- which(rownames(pairwise) == gene_symbol)[1]
  } else {
    stop("Could not identify gene names in pairwise object.")
  }

  if (is.na(gene_idx)) {
    stop(sprintf("Gene '%s' not found in pairwise results.", gene_symbol))
  }

  sig_cols <- grep(paste0("^", sig_prefix), colnames(pairwise), value = TRUE)

  res <- lapply(sig_cols, function(sig_col) {
    contrast <- sub(paste0("^", sig_prefix), "", sig_col)
    groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]

    if (length(groups) != 2) return(NULL)

    g1 <- groups[1]
    g2 <- groups[2]

    if (!(g1 %in% plot_celltypes && g2 %in% plot_celltypes)) return(NULL)

    pval <- as.numeric(pairwise[gene_idx, sig_col])

    tibble(
      group1 = g1,
      group2 = g2,
      p_value = pval,
      label = p_to_stars(pval)
    )
  }) %>%
    bind_rows()

  if (nrow(res) == 0) return(res)

  res %>%
    dplyr::filter(!is.na(label), p_value < alpha_thresh) %>%
    dplyr::mutate(
      xmin = match(group1, plot_celltypes),
      xmax = match(group2, plot_celltypes)
    ) %>%
    dplyr::arrange(pmin(xmin, xmax), pmax(xmin, xmax))
}

make_sig_df <- function(expr_df, sig_df) {
  if (is.null(sig_df) || nrow(sig_df) == 0) return(NULL)

  y_max <- max(expr_df$expr, na.rm = TRUE)
  y_range <- diff(range(expr_df$expr, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) y_range <- 1

  sig_df %>%
    dplyr::mutate(
      level = row_number(),
      y = y_max + level * 0.10 * y_range,
      y_text = y + 0.02 * y_range,
      y_tip = y - 0.02 * y_range
    )
}

plot_gene_violin <- function(gene_symbol) {
  expr_df <- get_gene_expr_df(sce, gene_symbol)
  sig_df <- get_sig_pairs_for_gene(pairwise, gene_symbol)
  sig_df <- make_sig_df(expr_df, sig_df)

  p <- ggplot(expr_df, aes(x = cell_type, y = expr, fill = cell_type)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      color = "black",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = celltype_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    labs(
      x = NULL,
      y = "logcounts",
      title = gene_symbol
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "italic")
    )

  if (!is.null(sig_df) && nrow(sig_df) > 0) {
    p <- p +
      geom_segment(
        data = sig_df,
        aes(x = xmin, xend = xmax, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = sig_df,
        aes(x = xmin, xend = xmin, y = y_tip, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = sig_df,
        aes(x = xmax, xend = xmax, y = y_tip, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_text(
        data = sig_df,
        aes(x = (xmin + xmax) / 2, y = y_text, label = label),
        inherit.aes = FALSE,
        size = 3.5
      )
  }

  p
}

# =========================================================
# Generate and save plots
# =========================================================
plots <- lapply(genes_to_plot, plot_gene_violin)
names(plots) <- genes_to_plot

combined_plot <- wrap_plots(plots, nrow = 1) +
  plot_annotation(tag_levels = "A")

  # =========================================================
# Save combined plot
# =========================================================
ggsave(
  filename = file.path(outdir, "violin_DRD1_celltypes_combined_TSHZ1.pdf"),
  plot = combined_plot,
  width = 11,
  height = 3.5
)