library("SpatialExperiment")
library("HDF5Array")
library("lobstr")
library("here")
library("sessioninfo")
library("Matrix")
library("dplyr")
library("ggplot2")
library("tidyr")
library("scales")
library("purrr")

# =========================================================
# Load data
# =========================================================
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)

spatial_enr_dir <- here(
  "processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers",
  "01_precast", "pseudobulk_capture_area", "final_clusters"
)
load(file.path(spatial_enr_dir, "model_results_precast_clusters.Rdata"))

pairwise <- modeling_results$pairwise

# =========================================================
# Keep only plotted domains
# =========================================================
all_domain_levels <- c(
  "D1 islands",
  "Endothelial/Ependymal",
  "Excitatory",
  "Inhibitory",
  "MSN_1",
  "MSN_2",
  "MSN_3",
  "WM"
)

plot_domain_levels <- c(
  "MSN_1",
  "MSN_2",
  "MSN_3",
  "D1 islands"
)

spe$spatial_domains <- factor(
  as.character(spe$spatial_domains),
  levels = all_domain_levels
)

# Mapping between display names in SPE and pairwise result names
domain_to_pairwise <- c(
  "D1 islands" = "D1.islands",
  "Endothelial/Ependymal" = "Endothelial.Ependymal",
  "Excitatory" = "Excitatory",
  "Inhibitory" = "Inhibitory",
  "MSN_1" = "MSN.1",
  "MSN_2" = "MSN.2",
  "MSN_3" = "MSN.3",
  "WM" = "WM"
)
pairwise_to_domain <- setNames(names(domain_to_pairwise), domain_to_pairwise)

# =========================================================
# Plot settings
# =========================================================
genes_to_plot <- c(
  "CALB1",
  "CARTPT",
  "TAC1",
  "PDYN",
  "ADORA2A",
  "PENK", 
  "DRD1",
  "DRD2"
)

sig_prefix <- "fdr_"
alpha_thresh <- 0.05

outdir <- here("plots", "25_revisions", "08-violin_plots_markers")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Colors for the 4 plotted domains
plot_cols <- c(
  "MSN_1" = "#66A61E",
  "MSN_2" = "#1B9E77",
  "MSN_3" = "#7570B3",
  "D1 islands" = "#E7298A"
)

# Only show D1 islands vs MSN comparisons
target_pairs <- tibble::tribble(
  ~group1,      ~group2,
  "MSN_1",      "D1 islands",
  "MSN_2",      "D1 islands",
  "MSN_3",      "D1 islands",
  "MSN_1",      "MSN_2",
  "MSN_1",      "MSN_3", 
  "MSN_2",       "MSN_3"
)

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

get_gene_expr_df <- function(spe, gene_symbol) {
  expr <- logcounts(spe)

  if (!gene_symbol %in% rownames(expr)) {
    stop(sprintf("Gene '%s' not found in rownames(logcounts(spe)).", gene_symbol))
  }

  data.frame(
    spatial_domain = as.character(spe$spatial_domains),
    expr = as.numeric(expr[gene_symbol, ]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(spatial_domain %in% plot_domain_levels) %>%
    dplyr::mutate(
      spatial_domain = factor(spatial_domain, levels = plot_domain_levels)
    )
}

get_sig_pairs_for_gene <- function(pairwise, gene_symbol) {
  gene_idx <- which(pairwise$gene == gene_symbol)[1]

  if (is.na(gene_idx)) {
    stop(sprintf("Gene '%s' not found in pairwise$gene.", gene_symbol))
  }

  sig_cols <- grep(paste0("^", sig_prefix), colnames(pairwise), value = TRUE)

  res <- lapply(sig_cols, function(sig_col) {
    contrast <- sub(paste0("^", sig_prefix), "", sig_col)
    groups <- strsplit(contrast, "-", fixed = TRUE)[[1]]

    if (length(groups) != 2) return(NULL)

    g1 <- pairwise_to_domain[[groups[1]]]
    g2 <- pairwise_to_domain[[groups[2]]]

    if (is.null(g1) || is.null(g2)) return(NULL)
    if (!(g1 %in% plot_domain_levels && g2 %in% plot_domain_levels)) return(NULL)

    pval <- as.numeric(pairwise[gene_idx, sig_col])

    data.frame(
      group1 = g1,
      group2 = g2,
      p_value = pval,
      label = p_to_stars(pval),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  if (nrow(res) == 0) return(res)

  res %>%
    dplyr::filter(!is.na(label), p_value < alpha_thresh) %>%
    dplyr::filter(
      (group1 == "D1 islands" & group2 %in% c("MSN_1", "MSN_2", "MSN_3")) |
      (group2 == "D1 islands" & group1 %in% c("MSN_1", "MSN_2", "MSN_3"))
    ) %>%
    dplyr::mutate(
      xmin = match(group1, plot_domain_levels),
      xmax = match(group2, plot_domain_levels)
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
  expr_df <- get_gene_expr_df(spe, gene_symbol)
  sig_df <- get_sig_pairs_for_gene(pairwise, gene_symbol)
  sig_df <- make_sig_df(expr_df, sig_df)

  p <- ggplot(expr_df, aes(x = spatial_domain, y = expr, fill = spatial_domain)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      color = "black",
      linewidth = 0.3
    ) +
    scale_fill_manual(values = plot_cols, drop = FALSE) +
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

for (g in genes_to_plot) {
  ggsave(
    filename = file.path(outdir, paste0("violin_", g, ".pdf")),
    plot = plots[[g]],
    width = 4.5,
    height = 3.0
  )
}

