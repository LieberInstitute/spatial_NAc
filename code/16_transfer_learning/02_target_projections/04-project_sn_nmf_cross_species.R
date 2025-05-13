rm(list = ls())
library(SingleCellExperiment)
library(HDF5Array)
library(scater)
library(scran)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcppML)
library(pheatmap)
library(ggrastr)
library(Seurat)
library(here)

# Load the rat case control data
opt <- list()
opt$data <- "rat_case_control_morphine_acute"
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
nmf_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", "human_NAc")

sce <- readRDS(file.path(dat_dir, "snRNA_seq_NAc.rds"))
orthologs_df <- readRDS(file.path(dat_dir, "orthologs_df.rds"))

# Read in the NMF
nmf_results <- readRDS(file.path(nmf_dir, "nmf_results.rds"))
loadings <- nmf_results@w

# Get the gene data from spe to map the ensembl ID to symbol
raw_in_path <- here("processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5")
spe <- loadHDF5SummarizedExperiment(raw_in_path)
geneData <- rowData(spe)
geneData <- geneData[geneData$gene_id %in% rownames(loadings), ]
geneData <- geneData[match(rownames(loadings), geneData$gene_id), ]
rownames(loadings) <- make.names(geneData$gene_name, unique = TRUE)

# Obtain the log counts from the rat data
logcounts <- sce[["RNA"]]$data

# Common genes
orthologs_df <- orthologs_df[orthologs_df$rat_genes %in% rownames(logcounts), ]
orthologs_df <- orthologs_df[orthologs_df$human_genes %in% rownames(loadings), ]
loadings <- loadings[rownames(loadings) %in% orthologs_df$human_genes, ]
logcounts <- logcounts[rownames(logcounts) %in% orthologs_df$rat_genes, ]
loadings <- loadings[match(orthologs_df$human_genes, rownames(loadings)), ]
logcounts <- logcounts[match(orthologs_df$rat_genes, rownames(logcounts)), ]

# Compute and add projections
proj <- project(w=loadings, data=logcounts)
proj1 <- proj[rowSums(proj) == 0, ]
proj2 <- proj[rowSums(proj) != 0,]

proj2 <- apply(proj2,1,function(x){x/sum(x)})
proj1 <- t(proj1)

proj_final <- cbind(proj2, proj1)
proj_final <- proj_final[ ,match(rownames(proj), colnames(proj_final))]

sce@meta.data <- cbind(sce@meta.data, proj_final)

if(opt$data == "rat_case_control_acute"){
    res_file <- here::here("processed-data", "16_transfer_learning", "02_target_projections", "projection_of_human_factors", "rat_case_control_cocaine_acute.rds")
}else{
    if(opt$data == "rat_case_control_repeated"){
    res_file <- here::here("processed-data", "16_transfer_learning", "02_target_projections", "projection_of_human_factors", "rat_case_control_cocaine_repeated.rds")
}else{
    res_file <- here::here("processed-data", "16_transfer_learning", "02_target_projections", "projection_of_human_factors", paste0(opt$data, ".rds"))
}
}

saveRDS(sce, res_file)



