####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)

opt <- list()
opt$data <- "human_NAc"
res_dir <- here("processed-data", "16_transfer_learning", "02_target_projections", opt$data)
plot_dir <- here("plots", "16_transfer_learning", "02_target_projections", opt$data)

spe <- readRDS(file = file.path(res_dir, paste0("spe_NMF.rds")))

compute_visium_gene_correlations <- function(spe, nmf_factors, assay_name = "logcounts") {
  expr <- assay(spe, assay_name)
  coldata <- colData(spe)

  cor_mat <- sapply(nmf_factors, function(fac) {
    vec <- coldata[[fac]]
    apply(expr, 1, function(g) cor(g, vec, method = "pearson"))
  })

  rownames(cor_mat) <- rownames(expr)
  return(cor_mat)
}


plot_visium_gene_correlation_heatmap <- function(cor_mat, nmf_factors, top_n = 20, output_dir, file_name) {
  library(dplyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)

  # 1. Select top_n genes per factor
  top_genes_list <- lapply(nmf_factors, function(fac) {
    if (!fac %in% colnames(cor_mat)) {
      stop(paste("Factor", fac, "not found in cor_mat"))
    }
    sort_genes <- sort(cor_mat[, fac], decreasing = TRUE)
    names(head(sort_genes, top_n))
  })
  names(top_genes_list) <- nmf_factors

  all_top_genes <- unique(unlist(top_genes_list))
  cor_mat_sub <- cor_mat[all_top_genes, nmf_factors, drop = FALSE]

  # Transpose to make rows = factors, columns = genes
  cor_mat_sub <- t(cor_mat_sub)

  # 2. Build gene-to-factor map (reverse roles for transposed matrix)
  gene_to_factor <- do.call(rbind, lapply(names(top_genes_list), function(fac) {
    genes <- top_genes_list[[fac]]
    if (length(genes) == 0) return(NULL)
    data.frame(factor = fac, gene = genes)
  }))

  if (is.null(gene_to_factor) || nrow(gene_to_factor) == 0) {
    stop("No top genes found for any NMF factor.")
  }

  gene_factor_labels <- gene_to_factor %>%
    group_by(gene) %>%
    summarise(Factor = paste(sort(unique(factor)), collapse = ", ")) %>%
    column_to_rownames("gene")

  # 3. Custom order of factor combinations
  all_factors <- sort(nmf_factors)
  all_combinations <- unlist(lapply(length(all_factors):1, function(i) {
    combn(all_factors, i, simplify = FALSE)
  }), recursive = FALSE)
  custom_levels <- sapply(all_combinations, function(facs) paste(facs, collapse = ", "))

  gene_factor_labels$Factor <- factor(gene_factor_labels$Factor, levels = custom_levels)
  gene_factor_labels <- gene_factor_labels[!is.na(gene_factor_labels$Factor), , drop = FALSE]
  cor_mat_sub <- cor_mat_sub[, rownames(gene_factor_labels), drop = FALSE]

  # 4. Define color scale and annotation
  n_groups <- length(levels(gene_factor_labels$Factor))
  palette_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_groups)
  factor_colors <- setNames(palette_colors, levels(gene_factor_labels$Factor))

  col_ha <- HeatmapAnnotation(
    Factor = gene_factor_labels$Factor,
    col = list(Factor = factor_colors),
    show_annotation_name = TRUE
  )

  # 5. Color mapping
  cor_range <- range(cor_mat_sub, na.rm = TRUE)
  breaks <- pretty(cor_range, n = 7)
  cor_colors <- circlize::colorRamp2(breaks, colorRampPalette(c("blue", "white", "red"))(length(breaks)))

  # Calculate dynamic font size for column labels
  n_genes <- ncol(cor_mat_sub)
  label_fontsize <- max(5, min(10, 300 / n_genes))

  # 6. Heatmap with genes as columns
  ht <- Heatmap(
    cor_mat_sub,
    name = "Correlation",
    col = cor_colors,
    top_annotation = col_ha,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = gene_factor_labels$Factor,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = label_fontsize),
    row_title = "NMF Factors",
    column_title = "Genes",
    heatmap_legend_param = list(title = "Correlation")
  )

  # 7. Save to file
  pdf(file.path(output_dir, file_name), width = 12, height = 4)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = "Gene Correlation with NMF Projections",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
}

plot_upset_top_correlated_genes <- function(cor_mat, nmf_factors, top_n = 200, output_dir, file_name = "upset_top_correlated_genes.pdf") {
  library(ComplexUpset)
  library(ggplot2)

  # 1. Get top N genes per factor
  top_gene_list <- lapply(nmf_factors, function(fac) {
    cor_values <- cor_mat[, fac]
    names(sort(cor_values, decreasing = TRUE)[1:top_n])
  })
  names(top_gene_list) <- nmf_factors

  # 2. Convert to binary membership matrix
  gene_df <- reshape2::melt(top_gene_list)
  colnames(gene_df) <- c("gene", "factor")
  gene_df$presence <- 1
  gene_matrix <- reshape2::acast(gene_df, gene ~ factor, fill = 0)
  gene_df_upset <- as.data.frame(gene_matrix)
  gene_df_upset$gene <- rownames(gene_df_upset)

  # 3. Plot with ComplexUpset
  pdf(file.path(output_dir, file_name), width = 9, height = 6)

  print(
    ComplexUpset::upset(
      gene_df_upset,
      intersect = nmf_factors,
      name = "Top 200 by correlation",
      sort_sets = "descending",
      matrix = (
        intersection_matrix(geom = geom_point(shape = 'circle filled', size = 3)) +
          scale_color_manual(
            values = c(nmf34 = "#D81B60", nmf35 = "#1E88E5", nmf44 = "#FFC107"),
            guide = guide_legend(override.aes = list(shape = 'circle'))
          )
      ),
      queries = list(
        upset_query(set = "nmf34", fill = "#D81B60"),
        upset_query(set = "nmf35", fill = "#1E88E5"),
        upset_query(set = "nmf44", fill = "#FFC107")
      )
    )
  )

  dev.off()
}


nmf_factors <- c("nmf34", "nmf35", "nmf44")
cor_mat <- compute_visium_gene_correlations(spe, nmf_factors)
rownames(cor_mat) <- make.unique(rowData(spe)$gene_name)

plot_visium_gene_correlation_heatmap(
  cor_mat = cor_mat,
  nmf_factors = nmf_factors,
  top_n = 20,
  output_dir = plot_dir,
  file_name = "visium_gene_projection_heatmap_D1_islands.pdf"
)

plot_upset_top_correlated_genes(
  cor_mat = cor_mat,
  nmf_factors = nmf_factors,
  top_n = 200,
  output_dir = plot_dir,
  file_name = "upset_top_200_correlated_genes.pdf"
)