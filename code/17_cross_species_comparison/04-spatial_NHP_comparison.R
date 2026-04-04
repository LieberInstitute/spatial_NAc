library(SpatialExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(Seurat)
library(RColorBrewer)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(getopt)
library(babelgene)
library(pheatmap)
library(ComplexHeatmap)

spe_dir <- here("processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
res_dir <- here("processed-data", "17_cross_species_comparison", "non_human_primate_spatial")
plot_dir <- here("plots", "17_cross_species_comparison", "non_human_primate_spatial")

subset_neurons <- FALSE
print(subset_neurons)

# -----------------------------
# Read in human 10x Visium
# -----------------------------
spe <- loadHDF5SummarizedExperiment(spe_dir)

# -----------------------------
# Read in NHP reference data
# -----------------------------
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")

if (subset_neurons) {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
} else {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}

# Change the hyphen to a period
if (subset_neurons) {
  Idents(sce_ext) <- sce_ext$MSN_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
} else {
  Idents(sce_ext) <- sce_ext$cell_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
}

# -----------------------------
# Load human spatial pseudobulk markers
# -----------------------------
reg_dir <- here(
  "processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers",
  "01_precast", "pseudobulk_capture_area", "final_clusters"
)

load(file.path(reg_dir, "model_results_precast_clusters.Rdata"))
sp_enrichment <- modeling_results$enrichment
rm(modeling_results)

if (subset_neurons) {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands")
} else {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory", "WM", "Endothelial.Ependymal")
}

# -----------------------------
# Add PRECAST labels to spe
# -----------------------------
clusters_file <- here(
  "processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast",
  "final_clusters", "precast_clusters.csv"
)

spe[["precast_clusters"]] <- colData(spe) |>
  as_tibble() |>
  left_join(read.csv(clusters_file), by = "key") |>
  pull(cluster) |>
  as.factor()

spe <- spe[, !is.na(spe[["precast_clusters"]])]

spe[["precast_clusters"]] <- gsub(" ", ".", spe[["precast_clusters"]])
spe[["precast_clusters"]] <- gsub("/", ".", spe[["precast_clusters"]])

spe <- spe[, spe$precast_clusters %in% select_domains]

# -----------------------------
# Build DEG table
# -----------------------------
DEG_Lists <- list()

for (i in select_domains) {
  df <- sp_enrichment[, grep(i, colnames(sp_enrichment)), drop = FALSE]
  colnames(df) <- gsub(paste0("_", i), "", colnames(df))
  df <- df[df$fdr < 0.05, , drop = FALSE]
  df <- df[df$logFC > 0, , drop = FALSE]
  df <- df[order(df$logFC, decreasing = TRUE), , drop = FALSE]
  df <- df[1:min(250, dim(df)[1]), , drop = FALSE]
  df$CellType <- i
  df$Gene_ID <- rownames(df)
  DEG_Lists[[i]] <- df
}

DEGs_df <- do.call(rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

geneData <- rowData(spe)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

# -----------------------------
# Map human -> macaque orthologs
# -----------------------------
Monkey_Orthos <- orthologs(
  genes = unique(DEGs_df$Gene_symbol),
  species = "Macaca mulatta"
)
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

DEGs_df <- merge(
  x = DEGs_df,
  y = Monkey_Orthos,
  by.x = "Gene_symbol",
  by.y = "human_symbol"
)

DEGs_df <- DEGs_df[!is.na(DEGs_df$Monkey_Symbol), , drop = FALSE]

# -----------------------------
# True 1:1 ortholog filtering at gene-pair level
# -----------------------------
ortholog_pairs <- unique(DEGs_df[, c("Gene_symbol", "Monkey_Symbol")])

human_counts <- table(ortholog_pairs$Gene_symbol)
nhp_counts <- table(ortholog_pairs$Monkey_Symbol)

valid_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% names(human_counts[human_counts == 1]) &
    ortholog_pairs$Monkey_Symbol %in% names(nhp_counts[nhp_counts == 1]),
  ,
  drop = FALSE
]

DEGs_df <- DEGs_df[
  DEGs_df$Gene_symbol %in% valid_pairs$Gene_symbol &
    DEGs_df$Monkey_Symbol %in% valid_pairs$Monkey_Symbol,
  ,
  drop = FALSE
]

# -----------------------------
# Restrict to genes present in NHP dataset
# -----------------------------
DEGs_df <- DEGs_df[
  DEGs_df$Monkey_Symbol %in% rownames(sce_ext[["RNA"]]$data),
  ,
  drop = FALSE
]

rownames(spe) <- rowData(spe)$gene_name

human_genes_use <- unique(DEGs_df$Gene_symbol)
nhp_genes_use <- unique(DEGs_df$Monkey_Symbol)

# -----------------------------
# Human spatial t-stat matrix
# -----------------------------
t_stat_mat_human <- matrix(
  nrow = length(human_genes_use),
  ncol = length(select_domains)
)
rownames(t_stat_mat_human) <- human_genes_use
colnames(t_stat_mat_human) <- select_domains

for (i in rownames(t_stat_mat_human)) {
  print(i)
  for (l in colnames(t_stat_mat_human)) {
    All_Others <- select_domains[which(select_domains != l)]

    c1_vals <- logcounts(spe)[i, spe$precast_clusters == l]
    c2_vals <- logcounts(spe)[i, spe$precast_clusters %in% All_Others]

    x_c1 <- mean(c1_vals)
    x_c2 <- mean(c2_vals)

    sdev_c1 <- sd(x = c1_vals)
    sdec_c2 <- sd(x = c2_vals)

    t_stat_mat_human[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) *
      sqrt(sum(spe$precast_clusters == l))
  }
}

if (subset_neurons) {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons_updated.rds"))
} else {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_updated.rds"))
}

# -----------------------------
# NHP t-stat matrix
# -----------------------------
t_stat_mat_nhp <- matrix(
  nrow = length(nhp_genes_use),
  ncol = length(levels(sce_ext))
)
rownames(t_stat_mat_nhp) <- nhp_genes_use
colnames(t_stat_mat_nhp) <- levels(sce_ext)

for (i in rownames(t_stat_mat_nhp)) {
  print(i)
  for (l in colnames(t_stat_mat_nhp)) {
    All_Others <- levels(sce_ext)[which(levels(sce_ext) != l)]

    c1_vals <- GetAssayData(object = sce_ext, assay = "RNA", slot = "data")[i, WhichCells(sce_ext, idents = l)]
    c2_vals <- GetAssayData(object = sce_ext, assay = "RNA", slot = "data")[i, WhichCells(sce_ext, idents = All_Others)]

    x_c1 <- mean(c1_vals)
    x_c2 <- mean(c2_vals)

    sdev_c1 <- sd(x = c1_vals)
    sdec_c2 <- sd(x = c2_vals)

    t_stat_mat_nhp[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) *
      sqrt(length(WhichCells(sce_ext, idents = l)))
  }
}

if (subset_neurons) {
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp_neurons_updated.rds"))
} else {
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp_updated.rds"))
}

# -----------------------------
# Correlation heatmap
# -----------------------------
t_stat_mat_nhp <- t_stat_mat_nhp[rowSums(is.na(t_stat_mat_nhp)) == 0, , drop = FALSE]
t_stat_mat_human <- t_stat_mat_human[rowSums(is.na(t_stat_mat_human)) == 0, , drop = FALSE]

ortholog_pairs <- unique(DEGs_df[, c("Gene_symbol", "Monkey_Symbol")])
ortholog_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% rownames(t_stat_mat_human) &
    ortholog_pairs$Monkey_Symbol %in% rownames(t_stat_mat_nhp),
  ,
  drop = FALSE
]

t_stat_mat_human <- t_stat_mat_human[
  match(ortholog_pairs$Gene_symbol, rownames(t_stat_mat_human)),
  ,
  drop = FALSE
]

t_stat_mat_nhp <- t_stat_mat_nhp[
  match(ortholog_pairs$Monkey_Symbol, rownames(t_stat_mat_nhp)),
  ,
  drop = FALSE
]

corr_mat <- matrix(NA, nrow = ncol(t_stat_mat_nhp), ncol = ncol(t_stat_mat_human))
rownames(corr_mat) <- colnames(t_stat_mat_nhp)
colnames(corr_mat) <- colnames(t_stat_mat_human)

for (i in rownames(corr_mat)) {
  for (j in colnames(corr_mat)) {
    corr_mat[i, j] <- cor(as.numeric(t_stat_mat_nhp[, i]), as.numeric(t_stat_mat_human[, j]))
  }
}

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("\\.", " ", colnames(corr_mat))

if (subset_neurons) {
  row_order <- c("D1 Matrix", "D1 Striosome", "D2 Matrix", "D2 Striosome", "D2 Shell/OT", "D1 Shell/OT", "D1/D2 Hybrid", "D1 ICj", "D1 NUDAP")
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands")
} else {
  rownames(corr_mat) <- gsub("_", " ", rownames(corr_mat))
  row_order <- c("DRD1", "DRD2", "Interneurons", "Astrocytes", "Microglia", "Oligos Pre", "Oligos", "Mural/Fibroblast", "Endothelial")
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial Ependymal")
}

corr_mat <- corr_mat[, match(col_order, colnames(corr_mat)), drop = FALSE]
corr_mat <- corr_mat[match(row_order, rownames(corr_mat)), , drop = FALSE]

if (subset_neurons) {
  saveRDS(corr_mat, file.path(res_dir, "corr_mat_nhp_neurons_spatial.rds"))
} else {
  saveRDS(corr_mat, file.path(res_dir, "corr_mat_nhp_spatial.rds"))
}

col_fun <- circlize::colorRamp2(c(-0.6, 0, 0.6), hcl_palette = "Purple-Green")

complex_plot_cor <- ComplexHeatmap::Heatmap(
  matrix = t(corr_mat),
  name = "Correlation",
  column_title = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = NULL,
  rect_gp = grid::gpar(col = "white", lwd = 2),
  col = col_fun,
  border = TRUE,
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = grid::unit(6, "cm"),
    title_position = "topcenter",
    title_gp = grid::gpar(fontsize = 14),
    border = "black",
    at = c(-0.65, -0.2, 0.2, 0.65)
  ),
  row_names_side = "left"
)

if (subset_neurons) {
  pdf(file.path(plot_dir, "correlation_heatmap_sn_neurons_updated.pdf"), width = 5, height = 4)
  ComplexHeatmap::draw(complex_plot_cor, heatmap_legend_side = "bottom")
  dev.off()
} else {
  pdf(file.path(plot_dir, "correlation_heatmap_sn_updated.pdf"), width = 7, height = 5)
  ComplexHeatmap::draw(complex_plot_cor, heatmap_legend_side = "bottom")
  dev.off()
}