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
#opt$subset_neurons <- FALSE

if (is.null(opt$subset_neurons)) {
  opt$subset_neurons <- FALSE
}

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control")

# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[, !sce$CellType.Final == "Neuron_Ambig"]

# Read in the reference data set
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

# Processing the human registration results
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")

# Load snRNA-seq registration results
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
DEGs_df <- do.call(what = rbind, DEG_Lists)
rownames(DEGs_df) <- seq_len(nrow(DEGs_df))

geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA_character_
for (i in seq_len(nrow(DEGs_df))) {
  DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]][1]
}

# Remove rows with missing human symbols
DEGs_df <- DEGs_df[!is.na(DEGs_df$Gene_symbol), , drop = FALSE]

# Add in the rat symbols
rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
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

# True 1:1 ortholog filtering at the unique gene-pair level
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

# Only select those genes that we have expression in the rat data for
if (opt$subset_neurons) {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext_neurons), , drop = FALSE]
} else {
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext), , drop = FALSE]
}

rownames(sce) <- rowData(sce)$gene_name

human_genes_use <- unique(DEGs_df$Gene_symbol)
rat_genes_use <- unique(DEGs_df$Rat_Gene_Symbol)

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

    t_stat_mat_human[i, l] <- (x_c1 - x_c2) / sqrt((sdev_c1^2 + sdec_c2^2) / 2) * sqrt(sum(sce$CellType.Final == l))
  }
}

if (opt$subset_neurons) {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons.rds"))
} else {
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human.rds"))
}

# Compute the t-stats for the rodent data
if (opt$subset_neurons) {
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

  saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat_neurons.rds"))
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

  saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat.rds"))
}