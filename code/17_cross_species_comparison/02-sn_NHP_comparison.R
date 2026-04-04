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

if (is.null(opt$subset_neurons)) {
  opt$subset_neurons <- FALSE
}

print(opt$subset_neurons)

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "non_human_primate")
plot_dir <- here("plots", "17_cross_species_comparison", "non_human_primate")
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")

# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]

# Read in the reference data set
if (opt$subset_neurons) {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
} else {
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}

if (opt$subset_neurons) {
  pdf(file.path(plot_dir, "QC_metrics_NHP_neurons.pdf"), width = 5, height = 4)
  print(
    VlnPlot(sce_ext, features = "nFeature_RNA", group.by = "MSN_type", pt.size = 0) +
      ggtitle("Number of genes") + theme(legend.position = "none") + xlab("")
  )
  dev.off()
} else {
  pdf(file.path(plot_dir, "QC_metrics_NHP.pdf"), width = 7, height = 4)
  print(
    VlnPlot(sce_ext, features = "nCount.RNA", group.by = "cell_type", pt.size = 0) +
      ggtitle("UMI count") + theme(legend.position = "none") + xlab("")
  )
  print(
    VlnPlot(sce_ext, features = "nFeatures_RNA", group.by = "cell_type", pt.size = 0) +
      ggtitle("Number of genes") + theme(legend.position = "none") + xlab("")
  )
  dev.off()
}

# Change the hyphen to a period
if (opt$subset_neurons) {
  Idents(sce_ext) <- sce_ext$MSN_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
} else {
  Idents(sce_ext) <- sce_ext$cell_type
  Idents(sce_ext) <- gsub(Idents(sce_ext), pattern = "-", replacement = ".")
}

# Processing the human registration results
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]

if (opt$subset_neurons) {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B"
  )
} else {
  select_cell_types <- c(
    "DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A",
    "DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
    "Astrocyte_A", "Astrocyte_B", "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial"
  )
}

sce <- sce[, sce$CellType.Final %in% select_cell_types]
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

# Create a dataframe
DEGs_df <- do.call(rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

# Remove rows with missing human symbols
DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

# Add in the macaque symbols
Monkey_Orthos <- orthologs(genes = unique(DEGs_df$Gene_symbol), species = "Macaca mulatta")
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

# Merge DEGs with macaque orthologs
DEGs_df <- merge(
  x = DEGs_df,
  y = Monkey_Orthos,
  by.x = "Gene_symbol",
  by.y = "human_symbol"
)

# Keep rows with mapped macaque symbols
DEGs_df <- DEGs_df[!is.na(DEGs_df$Monkey_Symbol), , drop = FALSE]

# True 1:1 ortholog filtering at the unique gene-pair level
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

# Keep only genes present in macaque data
DEGs_df <- DEGs_df[
  DEGs_df$Monkey_Symbol %in% rownames(sce_ext[["RNA"]]$data),
  ,
  drop = FALSE
]

rownames(sce) <- rowData(sce)$gene_name

human_genes_use <- unique(DEGs_df$Gene_symbol)
nhp_genes_use <- unique(DEGs_df$Monkey_Symbol)

t_stat_mat_human <- matrix(
  nrow = length(human_genes_use),
  ncol = length(select_cell_types)
)
rownames(t_stat_mat_human) <- human_genes_use
colnames(t_stat_mat_human) <- select_cell_types

for (i in rownames(t_stat_mat_human)) {
  print(i)
  for (l in colnames(t_stat_mat_human)) {
    All_Others <- select_cell_types[which(select_cell_types != l)]

    c1_vals <- logcounts(sce)[i, sce$CellType.Final == l]
    c2_vals <- logcounts(sce)[i, sce$CellType.Final %in% All_Others]

    x_c1 <- mean(c1_vals)
    x_c2 <- mean(c2_vals)

    sdev_c1 <- sd(x = c1_vals)
    sdec_c2 <- sd(x = c2_vals)

    t_stat_mat_human[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) *
      sqrt(sum(sce$CellType.Final == l))
  }
}

if (opt$subset_neurons) {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons.rds"))
} else {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human.rds"))
}

# Compute the t-stats for the NHP data
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

if (opt$subset_neurons) {
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp_neurons.rds"))
} else {
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp.rds"))
}