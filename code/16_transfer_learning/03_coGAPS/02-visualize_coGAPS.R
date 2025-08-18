rm(list = ls())
library(Seurat)
library(Matrix)
library(CoGAPS)
library(here)
library(ggplot2)
library(pheatmap)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(scran)
library(scuttle)
library(SingleCellExperiment)

opt <- list()
opt$data <- "rat_case_control_repeated"

res_dir <- here::here("processed-data", "16_transfer_learning", "03_coGAPS")
res_dir <- paste0(res_dir, "/", opt$data)

# Read in the processed data
sce <- readRDS(file.path(res_dir, "snRNA_seq.rds"))
cogaps_result <- readRDS(file.path(res_dir, "cogaps_results.rds"))

# Extract results from cogaps
A_matrix <- getPatternMatrix(cogaps_result)
P_matrix <- getFeatureLoadings(cogaps_result)

for (i in 1:ncol(A_matrix)) {
  sce[[paste0("CoGAPS_Pattern", i)]] <- A_matrix[, i]
}

sce$CellType_Stim <- paste0(sce$Combo_CellType,"_", sce$Stim)

plot_dir <- here::here("plots", "16_transfer_learning", "03_coGAPS")
plot_dir <- paste0(plot_dir, "/", opt$data)
pdf(file.path(plot_dir, "coGAPS_VlnPlot.pdf"), width = 20, height = 7)
p_list <- list()
for(i in c(1:ncol(A_matrix))){
p_list[[i]] <- VlnPlot(sce, features = paste0("CoGAPS_Pattern", i), group.by = "CellType_Stim", pt.size = 0)
print(p_list[[i]] + theme(legend.position = "none"))
}
dev.off()
