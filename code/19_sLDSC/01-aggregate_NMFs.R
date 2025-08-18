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
library(edgeR)
library(scran)
library(scuttle)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", "human_NAc")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", "human_NAc")

sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))
dat <- sce[["RNA"]]$data
x <- readRDS(file.path(res_dir,paste0("nmf_results.rds")))

# Read in the gene data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)
geneData <- rowData(spe)
rm(spe)
geneData <- geneData[geneData$gene_id %in% rownames(dat), ]
geneData <- geneData[match(rownames(dat), geneData$gene_id), ]
geneData <- geneData[!duplicated(geneData$gene_name), ]
geneData <- geneData[geneData$gene_type == "protein_coding", ]
geneData <- geneData[!grepl("^MT-", geneData$gene_name), ]
geneData <- geneData[!grepl("^RP[SL]", geneData$gene_name), ]

# Compute correlations
loadings <- x@h
dat <- as.matrix(dat)
corr_mat <- cor(t(dat),t(loadings))

corr_mat <- corr_mat[rownames(corr_mat) %in% geneData$gene_id, ]
geneData <- geneData[match(rownames(corr_mat), geneData$gene_id), ]
rownames(corr_mat) <- geneData$gene_name

res_dir <- here("processed-data", "19_sLDSC", "human_NAc_NMF", "input_files")
write.table(corr_mat, file.path(res_dir, 'NMF_correlations.tsv'), col.names = TRUE,
            row.names = TRUE, sep = "\t")