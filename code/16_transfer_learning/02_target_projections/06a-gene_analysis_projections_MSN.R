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

nmf_factors <- c("nmf3", "nmf4", "nmf7", "nmf10", "nmf39")
cor_mat <- compute_visium_gene_correlations(spe, nmf_factors)
rownames(cor_mat) <- make.unique(rowData(spe)$gene_name)

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
  cor_mat_sub <- t(cor_mat_sub)  # Rows = factors, Columns = genes

  # 2. Build gene-to-factor mapping
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

  # Keep only observed combinations (fix: no unused combinations)
  observed_levels <- c("nmf10, nmf4, nmf7", "nmf10, nmf3", "nmf10, nmf4", "nmf10, nmf7", 
  "nmf4, nmf7",  "nmf3, nmf4", "nmf10", "nmf4", "nmf7", "nmf39", "nmf3")
  gene_factor_labels$Factor <- factor(gene_factor_labels$Factor, levels = observed_levels)
  gene_factor_labels <- gene_factor_labels[!is.na(gene_factor_labels$Factor), , drop = FALSE]
  cor_mat_sub <- cor_mat_sub[, rownames(gene_factor_labels), drop = FALSE]

  # 4. Define annotation colors
  n_groups <- length(levels(gene_factor_labels$Factor))
  palette_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_groups)
  factor_colors <- setNames(palette_colors, levels(gene_factor_labels$Factor))

  col_ha <- HeatmapAnnotation(
    Factor = gene_factor_labels$Factor,
    col = list(Factor = factor_colors),
    show_annotation_name = TRUE
  )

  # 5. Color mapping for correlation values
  cor_range <- range(cor_mat_sub, na.rm = TRUE)
  breaks <- pretty(cor_range, n = 7)
  cor_colors <- circlize::colorRamp2(breaks, colorRampPalette(c("blue", "white", "red"))(length(breaks)))

  # 6. Adjust font size dynamically
  n_genes <- ncol(cor_mat_sub)
  label_fontsize <- 10

  # 7. Draw heatmap
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

  # 8. Export to PDF
  pdf(file.path(output_dir, file_name), width = 14, height = 4)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       column_title = "Gene Correlation with NMF Projections",
       column_title_gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
}

plot_visium_gene_correlation_heatmap(
  cor_mat = cor_mat,
  nmf_factors = nmf_factors,
  top_n = 20,
  output_dir = plot_dir,
  file_name = "visium_gene_projection_heatmap_MSN.pdf"
)


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
  pdf(file.path(output_dir, file_name), width = 8, height = 4)

  print(
    ComplexUpset::upset(
      gene_df_upset,
      intersect = nmf_factors,
      name = "Top 200 by correlation",
      sort_sets = "descending",
      matrix = (
        intersection_matrix(geom = geom_point(shape = 'circle filled', size = 3)) +
          scale_color_manual(
            values = c(nmf3 = "#0072B2", nmf4 = "#D55E00", nmf7 = "#F0E442", nmf10 = "#009E73" , nmf39 = "#CC79A7"),
            guide = guide_legend(override.aes = list(shape = 'circle'))
          )
      ),
      base_annotations=list(
        'Intersection size'=intersection_size(counts=FALSE)
      ),
      queries = list(
        upset_query(set = "nmf3", fill = "#0072B2"),
        upset_query(set = "nmf4", fill = "#D55E00"),
        upset_query(set = "nmf7", fill = "#F0E442"), 
        upset_query(set = "nmf10", fill = "#009E73"), 
        upset_query(set = "nmf39", fill = "#CC79A7")
      )
    )
  )

  dev.off()
}


plot_upset_top_correlated_genes(
  cor_mat = cor_mat,
  nmf_factors = nmf_factors,
  top_n = 200,
  output_dir = plot_dir,
  file_name = "upset_top_200_correlated_genes_MSN.pdf"
)