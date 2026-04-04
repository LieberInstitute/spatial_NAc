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
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control_spatial")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control_spatial")

subset_neurons <- FALSE
t_stat_human_processed <- TRUE
t_stat_rat_processed <- FALSE
print(subset_neurons)

# -----------------------------
# Read in human 10x Visium
# -----------------------------
spe <- loadHDF5SummarizedExperiment(spe_dir)

# -----------------------------
# Read in external data
# -----------------------------
dat_dir <- here("processed-data", "12_snRNA")
sce_ext <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
sce_ext_neurons <- subset(
  sce_ext,
  idents = c(
    "Drd1-MSN-1", "Drd1-MSN-2", "Drd2-MSN-1", "Drd2-MSN-2",
    "Drd3-MSN", "Grm8-MSN", "GABAergic", "Chat-Interneuron",
    "Pvalb-Interneuron", "Sst-Interneuron", "Glutamatergic"
  )
)

# Change the hyphen to a period
Idents(sce_ext_neurons) <- gsub(Idents(sce_ext_neurons), pattern = "-", replacement = ".")
Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")

# -----------------------------
# Processing the human spatial registration results
# -----------------------------
reg_dir <- here(
  "processed-data", "10_post_clustering_analysis", "01_pseudobulk_markers",
  "01_precast", "pseudobulk_capture_area", "final_clusters"
)

load(file.path(reg_dir, "model_results_precast_clusters.Rdata"))
sp_enrichment <- modeling_results$enrichment
rm(modeling_results)

if (subset_neurons) {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory")
} else {
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", "Inhibitory", "WM", "Endothelial.Ependymal")
}

# -----------------------------
# Add PRECAST results to the spe object
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

# Remove spots with no PRECAST output
spe <- spe[, !is.na(spe[["precast_clusters"]])]

# Replace any spatial characters with a "."
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

DEGs_df <- do.call(what = rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

geneData <- rowData(spe)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

# Remove rows with missing human symbols
DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

# -----------------------------
# Add in rat ortholog symbols
# -----------------------------
rat_human_ortholog <- read.delim(
  "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt"
)
hom_rat <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
hom_hs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]

DEGs_df$Rat_Gene_Symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  print(i)
  human_keys <- hom_hs[hom_hs$Symbol %in% DEGs_df[i, "Gene_symbol"], "DB.Class.Key"]
  rat_symbols <- unique(hom_rat[hom_rat$DB.Class.Key %in% human_keys, "Symbol"])

  if (length(rat_symbols) == 1) {
    DEGs_df[i, "Rat_Gene_Symbol"] <- rat_symbols
  }
}

# Keep rows with mapped rat orthologs
DEGs_df <- DEGs_df[!is.na(DEGs_df$Rat_Gene_Symbol), , drop = FALSE]

# -----------------------------
# True 1:1 ortholog filtering at the unique gene-pair level
# -----------------------------
ortholog_pairs <- unique(DEGs_df[, c("Gene_symbol", "Rat_Gene_Symbol")])

human_counts <- table(ortholog_pairs$Gene_symbol)
rat_counts <- table(ortholog_pairs$Rat_Gene_Symbol)

valid_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% names(human_counts[human_counts == 1]) &
    ortholog_pairs$Rat_Gene_Symbol %in% names(rat_counts[rat_counts == 1]),
  ,
  drop = FALSE
]

DEGs_df <- DEGs_df[
  DEGs_df$Gene_symbol %in% valid_pairs$Gene_symbol &
    DEGs_df$Rat_Gene_Symbol %in% valid_pairs$Rat_Gene_Symbol,
  ,
  drop = FALSE
]

# -----------------------------
# Keep only genes present in the rat data
# -----------------------------
if (subset_neurons) {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext_neurons), , drop = FALSE]
} else {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext), , drop = FALSE]
}

rownames(spe) <- rowData(spe)$gene_name

human_genes_use <- unique(DEGs_df$Gene_symbol)
rat_genes_use <- unique(DEGs_df$Rat_Gene_Symbol)

# -----------------------------
# Human spatial t-stat matrix
# -----------------------------
if(t_stat_human_processed){
  if (subset_neurons) {
  t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human_neurons_updated.rds"))
  } else {
  t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human_updated.rds"))
  }

}else{
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

      t_stat_mat_human[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) * sqrt(sum(spe$precast_clusters == l))
    }
  }

  if (subset_neurons) {
    saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons_updated.rds"))
  } else {
    saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_updated.rds"))
  }
}



# -----------------------------
# Rodent t-stat matrix
# -----------------------------
if(t_stat_rat_processed){
  if (subset_neurons) {
  t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat_neurons_updated.rds"))
  } else {
  t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat_updated.rds"))
  }
}else{
    if (subset_neurons) {
    t_stat_mat_rat <- matrix(
      nrow = length(rat_genes_use),
      ncol = length(levels(sce_ext_neurons))
    )
    rownames(t_stat_mat_rat) <- rat_genes_use
    colnames(t_stat_mat_rat) <- levels(sce_ext_neurons)

    for (i in rownames(t_stat_mat_rat)) {
      print(i)
      for (l in colnames(t_stat_mat_rat)) {
        All_Others <- levels(sce_ext_neurons)[which(levels(sce_ext_neurons) != l)]

        c1_vals <- GetAssayData(object = sce_ext_neurons, assay = "RNA", slot = "data")[i, WhichCells(sce_ext_neurons, idents = l)]
        c2_vals <- GetAssayData(object = sce_ext_neurons, assay = "RNA", slot = "data")[i, WhichCells(sce_ext_neurons, idents = All_Others)]

        x_c1 <- mean(c1_vals)
        x_c2 <- mean(c2_vals)

        sdev_c1 <- sd(x = c1_vals)
        sdec_c2 <- sd(x = c2_vals)

        t_stat_mat_rat[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) * sqrt(length(WhichCells(sce_ext_neurons, idents = l)))
      }
    }

    saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat_neurons_updated.rds"))
  } else {
    t_stat_mat_rat <- matrix(
      nrow = length(rat_genes_use),
      ncol = length(levels(sce_ext))
    )
    rownames(t_stat_mat_rat) <- rat_genes_use
    colnames(t_stat_mat_rat) <- levels(sce_ext)

    for (i in rownames(t_stat_mat_rat)) {
      print(i)
      for (l in colnames(t_stat_mat_rat)) {
        All_Others <- levels(sce_ext)[which(levels(sce_ext) != l)]

        c1_vals <- GetAssayData(object = sce_ext, assay = "RNA", slot = "data")[i, WhichCells(sce_ext, idents = l)]
        c2_vals <- GetAssayData(object = sce_ext, assay = "RNA", slot = "data")[i, WhichCells(sce_ext, idents = All_Others)]

        x_c1 <- mean(c1_vals)
        x_c2 <- mean(c2_vals)

        sdev_c1 <- sd(x = c1_vals)
        sdec_c2 <- sd(x = c2_vals)

        t_stat_mat_rat[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) * sqrt(length(WhichCells(sce_ext, idents = l)))
      }
    }

    saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat_updated.rds"))
  }

}

# -----------------------------
# Correlation heatmap
# -----------------------------
t_stat_mat_rat <- t_stat_mat_rat[rowSums(is.na(t_stat_mat_rat)) == 0, , drop = FALSE]
t_stat_mat_human <- t_stat_mat_human[rowSums(is.na(t_stat_mat_human)) == 0, , drop = FALSE]

ortholog_pairs <- unique(DEGs_df[, c("Gene_symbol", "Rat_Gene_Symbol")])
ortholog_pairs <- ortholog_pairs[
  ortholog_pairs$Gene_symbol %in% rownames(t_stat_mat_human) &
    ortholog_pairs$Rat_Gene_Symbol %in% rownames(t_stat_mat_rat),
  ,
  drop = FALSE
]

t_stat_mat_human <- t_stat_mat_human[
  match(ortholog_pairs$Gene_symbol, rownames(t_stat_mat_human)),
  ,
  drop = FALSE
]

t_stat_mat_rat <- t_stat_mat_rat[
  match(ortholog_pairs$Rat_Gene_Symbol, rownames(t_stat_mat_rat)),
  ,
  drop = FALSE
]

corr_mat <- matrix(NA, nrow = ncol(t_stat_mat_rat), ncol = ncol(t_stat_mat_human))
rownames(corr_mat) <- colnames(t_stat_mat_rat)
colnames(corr_mat) <- colnames(t_stat_mat_human)

for (i in rownames(corr_mat)) {
  for (j in colnames(corr_mat)) {
    corr_mat[i, j] <- cor(as.numeric(t_stat_mat_rat[, i]), as.numeric(t_stat_mat_human[, j]))
  }
}

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("\\.", " ", colnames(corr_mat))

if (subset_neurons) {
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory")
  row_order <- c(
    "Drd1 MSN 1", "Drd2 MSN 1", "Drd1 MSN 2", "Drd2 MSN 2", "Drd3 MSN",
    "Grm8 MSN", "Glutamatergic", "GABAergic", "Pvalb Interneuron",
    "Chat Interneuron", "Sst Interneuron"
  )
} else {
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial Ependymal")
  row_order <- c(
    "Drd1 MSN 1", "Drd2 MSN 1", "Drd1 MSN 2", "Drd2 MSN 2", "Drd3 MSN", "Grm8 MSN",
    "Glutamatergic", "GABAergic", "Pvalb Interneuron", "Chat Interneuron", "Sst Interneuron",
    "Astrocyte", "Microglia", "Polydendrocyte", "Olig 1", "Mural"
  )
}

corr_mat <- corr_mat[, match(col_order, colnames(corr_mat)), drop = FALSE]
corr_mat <- corr_mat[match(row_order, rownames(corr_mat)), , drop = FALSE]

col_fun <- circlize::colorRamp2(c(-0.65, 0, 0.65), hcl_palette = "Purple-Green")

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