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

spec <- matrix(
    c(
        "subset_neurons", "n", 1, "logical", "Perform correlation only for neuronal cell types?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
#opt <- list()
#opt$subset_neurons <- TRUE

if (is.null(opt$subset_neurons)) {
  opt$subset_neurons <- FALSE
}

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control")

# -----------------------------
# Read in human snRNA-seq
# -----------------------------
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]

# -----------------------------
# Read in the reference data set
# -----------------------------
sce_ext <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
sce_ext_neurons <- subset(
  sce_ext,
  idents = c(
    "Drd1-MSN-1", "Drd1-MSN-2", "Drd2-MSN-1", "Drd2-MSN-2",
    "Drd3-MSN", "Grm8-MSN", "GABAergic", "Chat-Interneuron",
    "Pvalb-Interneuron", "Sst-Interneuron", "Glutamatergic"
  )
)

Idents(sce_ext_neurons) <- gsub(Idents(sce_ext_neurons), pattern = "-", replacement = ".")
Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")

# -----------------------------
# Human registration results
# -----------------------------
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]

if (opt$subset_neurons) {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F"
  )
} else {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
    "Astrocyte_A", "Astrocyte_B", "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial"
  )
}

sce <- sce[, sce$CellType.Final %in% select_cell_types]

# -----------------------------
# Build DEG table
# -----------------------------
DEG_Lists <- list()

for (i in select_cell_types) {
  df <- sn_enrichment[, grep(i, colnames(sn_enrichment)), drop = FALSE]
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

geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

# -----------------------------
# Map human -> rat orthologs
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

DEGs_df <- DEGs_df[!is.na(DEGs_df$Rat_Gene_Symbol), , drop = FALSE]

# -----------------------------
# True 1:1 ortholog filtering at gene-pair level
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
# Restrict to genes present in rodent dataset
# -----------------------------
if (opt$subset_neurons) {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext_neurons), , drop = FALSE]
} else {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext), , drop = FALSE]
}

rownames(sce) <- rowData(sce)$gene_name

# -----------------------------
# Read updated t-stat matrices
# -----------------------------
if (opt$subset_neurons) {
  t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human_neurons.rds"))
  t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat_neurons.rds"))
} else {
  t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human.rds"))
  t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat.rds"))
}

# -----------------------------
# Remove rows with NA values
# -----------------------------
t_stat_mat_rat <- t_stat_mat_rat[rowSums(is.na(t_stat_mat_rat)) == 0, , drop = FALSE]
t_stat_mat_human <- t_stat_mat_human[rowSums(is.na(t_stat_mat_human)) == 0, , drop = FALSE]

# -----------------------------
# Align using unique ortholog pairs
# -----------------------------
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

# -----------------------------
# Correlation matrix
# -----------------------------
corr_mat <- matrix(NA, nrow = ncol(t_stat_mat_rat), ncol = ncol(t_stat_mat_human))
rownames(corr_mat) <- colnames(t_stat_mat_rat)
colnames(corr_mat) <- colnames(t_stat_mat_human)

for (i in rownames(corr_mat)) {
  for (j in colnames(corr_mat)) {
    corr_mat[i, j] <- cor(as.numeric(t_stat_mat_rat[, i]), as.numeric(t_stat_mat_human[, j]))
  }
}

# -----------------------------
# Plot heatmap
# -----------------------------
colrange <- seq(-0.6, 0.6, by = 0.01)
colorpal <- colorRampPalette(brewer.pal(n = 7, name = "PRGn"))(length(colrange))

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("_", " ", colnames(corr_mat))

if (opt$subset_neurons) {
  row_order <- c(
    "Drd1 MSN 1", "Drd1 MSN 2", "Drd2 MSN 1", "Drd2 MSN 2", "Grm8 MSN",
    "Glutamatergic", "Pvalb Interneuron", "Chat Interneuron",
    "Sst Interneuron", "GABAergic", "Drd3 MSN"
  )
  col_order <- c(
    "DRD1 MSN A", "DRD1 MSN C", "DRD2 MSN A", "DRD2 MSN B", "DRD1 MSN B",
    "DRD1 MSN D", "Excitatory", "Inh A", "Inh B", "Inh C", "Inh D", "Inh E", "Inh F"
  )
} else {
  row_order <- c(
    "Drd1 MSN 1", "Drd1 MSN 2", "Drd2 MSN 1", "Drd2 MSN 2", "Grm8 MSN",
    "Glutamatergic", "Pvalb Interneuron", "Chat Interneuron",
    "Sst Interneuron", "GABAergic", "Drd3 MSN",
    "Astrocyte", "Microglia", "Polydendrocyte", "Olig 1", "Mural"
  )
  col_order <- c(
    "DRD1 MSN A", "DRD1 MSN C", "DRD2 MSN A", "DRD2 MSN B", "DRD1 MSN B",
    "DRD1 MSN D", "Excitatory", "Inh A", "Inh B", "Inh C", "Inh D", "Inh E", "Inh F",
    "Astrocyte A", "Astrocyte B", "Microglia", "OPC", "Oligo", "Endothelial", "Ependymal"
  )
}

corr_mat <- corr_mat[, match(col_order, colnames(corr_mat)), drop = FALSE]
corr_mat <- corr_mat[match(row_order, rownames(corr_mat)), , drop = FALSE]

if (opt$subset_neurons) {
  pdf(file.path(plot_dir, "correlation_heatmap_sn_neurons_updated.pdf"), width = 8, height = 5)
  pheatmap::pheatmap(
    corr_mat,
    color = colorpal,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    breaks = colrange,
    fontsize = 11,
    fontsize_row = 11,
    fontsize_col = 11,
    fontsize_number = 8,
    legend_breaks = seq(-0.6, 0.6, by = 0.3),
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black"
  )
  dev.off()
} else {
  pdf(file.path(plot_dir, "correlation_heatmap_sn_updated.pdf"), width = 10, height = 6)
  pheatmap::pheatmap(
    corr_mat,
    color = colorpal,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    breaks = colrange,
    fontsize = 11,
    fontsize_row = 11,
    fontsize_col = 11,
    fontsize_number = 8,
    legend_breaks = seq(-0.6, 0.6, by = 0.3),
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black"
  )
  dev.off()
}