rm(list = ls())

# ==== LIBRARIES ====
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(RcppML)
library(SingleCellExperiment)
library(ggrepel)
library(here)
library(tidyr)
library(ComplexUpset)
library(hrbrthemes)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)

# ==== LOAD DATA ====
sc_dir <- here::here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(sc_dir, "sce_CellType_noresiduals.Rds"))
geneData <- rowData(sce)

dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", "human_NAc")
sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

nmf_dir <- here::here('processed-data', '16_transfer_learning')
nmf <- readRDS(file.path(nmf_dir, "01_process_reference", "RCppML", "human_NAc", paste0("nmf_results.rds")))
H <- nmf@h
W <- nmf@w
expr <- sce[["RNA"]]$data

df <- geneData[geneData$gene_id %in% rownames(expr), ]
df <- df[match(rownames(expr), df$gene_id), ]
rownames(expr) <- make.unique(df$gene_name)
rownames(W) <- make.unique(df$gene_name)

factors <- c("nmf3", "nmf4", "nmf7", "nmf10", "nmf39")
top_n <- 30
upset_n <- 200
ranking_method <- "correlation"  # "loading", "correlation", or "combined"

output_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "gene_analysis", "human_NAc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== FUNCTIONS ====

get_top_genes_by_method <- function(scores_df, top_n, method = "combined") {
  stopifnot(method %in% c("loading", "correlation", "combined"))

  if (method == "loading") {
    scores_df <- scores_df %>% arrange(desc(loading))
  } else if (method == "correlation") {
    scores_df <- scores_df %>% arrange(desc(correlation))
  } else {
    scores_df <- scores_df %>%
      mutate(rank_score = rank(-loading) + rank(-correlation)) %>%
      arrange(rank_score)
  }

  return(head(scores_df$gene, top_n))
}

plot_top_loadings <- function(w, factor, top_n, outdir) {
  library(ggplot2)
  library(dplyr)
  library(scales)

  df <- data.frame(
    gene = rownames(w),
    loading = as.numeric(w[, factor])
  ) %>%
    dplyr::arrange(desc(loading)) %>%
    dplyr::slice(1:top_n)

  df$gene <- factor(df$gene, levels = df$gene[order(df$loading)])

  p <- ggplot(df, aes(x = gene, y = loading, fill = loading)) +
    geom_bar(stat = "identity") +
    scale_fill_gradientn(colors = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF")) +
    coord_flip() +
    labs(
      title = paste("Top", top_n, "Loadings for", factor),
      x = NULL,
      y = "Loading"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      panel.grid.major.y = element_blank(),   # Remove horizontal grid lines
      panel.grid.minor.y = element_blank(),   # Remove minor horizontal grid lines
      panel.grid.major.x = element_line(color = "grey80")  # Keep vertical grid lines
    )

  ggsave(file.path(outdir, paste0("barplot_", factor, ".pdf")), p, width = 5, height = 5)
}

compute_gene_cor <- function(expr, h_vec) {
  apply(expr, 1, function(gene_expr) cor(gene_expr, h_vec, method = "pearson"))
}

plot_loading_vs_correlation <- function(w, cors, factor, top_n, outdir) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  df <- data.frame(
    gene = rownames(w),
    loading = as.numeric(w[, factor]),
    correlation = cors
  )

  # Rank-based combined score
  df <- df %>%
    mutate(
      rank_loading     = rank(-loading),
      rank_correlation = rank(-correlation),
      rank_sum         = rank_loading + rank_correlation
    )

  # Top gene identification
  df <- df %>%
    mutate(
      top_loading     = rank_loading <= top_n,
      top_correlation = rank_correlation <= top_n
    ) %>%
    arrange(rank_sum) %>%
    mutate(top_combined = row_number() <= top_n)

  # Highlight category
  df$highlight <- "none"
  df$highlight[df$top_loading]     <- "loading"
  df$highlight[df$top_correlation] <- "correlation"
  df$highlight[df$top_combined]    <- "combined"
  df$highlight <- factor(df$highlight, levels = c("none", "loading", "correlation", "combined"))

  # Colorblind-friendly palette
  highlight_colors <- c(
    "none"        = "#CCCCCC",  # Gray
    "loading"     = "#0072B2",  # Blue
    "correlation" = "#E69F00",  # Orange
    "combined"    = "#8E0152"   # Dark magenta ("Vermillion")
  )

  # Plot with labeled top_combined genes
  p <- ggplot(df, aes(x = loading, y = correlation, color = highlight)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_text_repel(
      data = df[df$top_combined, ],
      aes(label = gene),
      size = 3,
      box.padding = 0.4,
      max.overlaps = 100,
      color = "black"
    ) +
    scale_color_manual(values = highlight_colors) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Loading vs Correlation:", factor),
      x = "Loading (W)", y = "Correlation with H",
      color = "Top Genes"
    )

  ggsave(file.path(outdir, paste0("scatter_loading_vs_corr_", factor, ".pdf")), p, width = 6, height = 6)

  return(df)
}

plot_gene_factor_correlation_complex <- function(expr, H, top_gene_list, output_dir, 
                                                 file_name = "gene_factor_correlation_heatmap_balanced.pdf", 
                                                 max_genes = 100) {
  library(circlize)
  library(ComplexHeatmap)

  # 1. Filter top genes and remove unwanted
  top_genes <- unique(unlist(top_gene_list))
  top_genes <- setdiff(top_genes, exclude_genes)
  expr_sub <- expr[top_genes, , drop = FALSE]

  # 2. Compute correlation matrix (genes x factors)
  selected_factors <- names(top_gene_list)
  cor_mat <- sapply(selected_factors, function(fac) {
    h_vec <- H[fac, ]
    apply(expr_sub, 1, function(g) cor(g, h_vec, method = "pearson"))
  })
  rownames(cor_mat) <- rownames(expr_sub)

  # 3. Balanced gene selection: highlight core and outlier
  core_factors <- c("nmf3", "nmf4", "nmf7", "nmf10")
  outlier_factors <- c("nmf39")

  core_var <- apply(cor_mat[, core_factors, drop = FALSE], 1, var)
  core_vs_outlier_diff <- abs(
    rowMeans(cor_mat[, core_factors, drop = FALSE]) - rowMeans(cor_mat[, outlier_factors, drop = FALSE])
  )
  combined_score <- scale(core_var) + scale(core_vs_outlier_diff)

  top_genes <- names(sort(combined_score[, 1], decreasing = TRUE))[1:min(max_genes, length(combined_score))]
  cor_mat <- cor_mat[top_genes, , drop = FALSE]

  # 4. Cluster rows (genes) by similarity
  row_order <- hclust(dist(cor_mat))$order

  # 5. Manually set column (factor) order
  desired_col_order <- c("nmf39", "nmf10", "nmf4", "nmf3", "nmf7")
  col_order <- desired_col_order[desired_col_order %in% colnames(cor_mat)]
  cor_mat <- cor_mat[, col_order, drop = FALSE]

  # 6. Color scale
  heatmap_colors <- circlize::colorRamp2(
    breaks = c(-0.2, 0, 0.3, 0.6),
    colors = c("navy", "white", "#fdae61", "#d73027")
  )

  # 7. Draw heatmap
  ht <- Heatmap(
    matrix = cor_mat,
    name = "Correlation",
    col = heatmap_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = row_order,
    column_order = col_order,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_title = "Genes",
    column_title = "NMF Factors",
    border = TRUE,
    heatmap_legend_param = list(
      title = "Correlation",
      at = c(-0.2, 0, 0.3, 0.6),
      labels = c("-0.2", "0", "0.3", "0.6")
    )
  )

  # 8. Save to PDF
  pdf(file.path(output_dir, file_name), width = 7, height = 10)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}


plot_upset_top_genes <- function(all_scores_list, output_dir, upset_n,
                                 method = "combined", file_name = "upset_plot.pdf") {
  top_gene_list <- lapply(all_scores_list, function(df) {
    get_top_genes_by_method(df, upset_n, method)
  })

  gene_df <- reshape2::melt(top_gene_list)
  colnames(gene_df) <- c("gene", "factor")
  gene_df$presence <- 1
  gene_matrix <- reshape2::acast(gene_df, gene ~ factor, fill = 0)
  gene_df_upset <- as.data.frame(gene_matrix)

  if (nrow(gene_df_upset) > 0 && ncol(gene_df_upset) > 1) {
     upset_plot <-  ComplexUpset::upset(
      gene_df_upset,
      intersect = colnames(gene_df_upset),
      name = "NMF",
      sort_sets = "descending", # MUST be a string in v1.3.6
      height_ratio = 0.6, 
      matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=3))
        + scale_color_manual(
            values=c('nmf34'='#D81B60', 'nmf35'='#1E88E5', 'nmf44'='#FFC107'),
            guide=guide_legend(override.aes=list(shape='circle'))
        )
      ),
      queries = list(upset_query(set = "nmf34", fill = "#D81B60"), 
                     upset_query(set = "nmf35", fill = "#1E88E5"), 
                     upset_query(set = "nmf44", fill = "#FFC107")))
    

    pdf(file.path(output_dir, file_name), width = 9, height = 6)
    print(upset_plot)
    dev.off()
  } else {
    warning("Insufficient data for UpSet plot.")
  }
}

# ==== MAIN ANALYSIS ====
top_combined_list <- list()
all_scores_list <- list()

for (fac in factors) {
  message("Processing ", fac)

  plot_top_loadings(W, fac, top_n, output_dir)

  h_vec <- H[fac, ]
  gene_cors <- compute_gene_cor(expr, h_vec)

  summary_df <- plot_loading_vs_correlation(W, gene_cors, fac, top_n, output_dir)

  write.csv(summary_df,
            file.path(output_dir, paste0("all_gene_scores_", fac, ".csv")),
            row.names = FALSE)

  selected_genes <- get_top_genes_by_method(summary_df, top_n, method = ranking_method)

  top_combined_list[[fac]] <- selected_genes
  all_scores_list[[fac]] <- summary_df
}

plot_gene_factor_correlation_complex(
  expr = expr,
  H = H,
  top_gene_list = top_combined_list,
  output_dir = output_dir,
  file_name = paste0("faceted_heatmap_", ranking_method, "_MSNs.pdf"),
  max_genes = 50
)

plot_upset_top_genes(
  all_scores_list,
  output_dir,
  upset_n,
  method = ranking_method,
  file_name = paste0("upset_top_genes_", ranking_method, "_MSNs.pdf")
)

message("All outputs saved to: ", output_dir)


