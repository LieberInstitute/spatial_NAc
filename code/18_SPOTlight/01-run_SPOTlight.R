library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(here)
library(HDF5Array)
library(dplyr)
library(tidyverse)
library(getopt)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(Matrix)
library(spatialNAcUtils)
library(spatialLIBD)

# Read in the annotated single nucleus data
sn_dir <- here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(sn_dir, "sce_CellType_noresiduals.Rds"))
sce$CellType.Final[sce$CellType.Final == "T-Cell"] <- "T_cell"

# Read in spatial data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Specify the output directory to save files to
out_path <- here("processed-data", "18_SPOTlight")

# Rename rows to symbols and make unique
rownames(sce) <- rowData(sce)$gene_name
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
rownames(sce) <- make.names(rownames(sce), unique = TRUE)

# Processing single nucleus data to obtain NMF factors
sce <- logNormCounts(sce)
# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^RP[SL]|MT\\.", x = rownames(sce))

dec <- modelGeneVar(sce, subset.row = genes)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
# Get the top 3000 HVGs
# Get the top 3000 genes.
hvg <- getTopHVGs(dec, n = 3000)

# Obtain marker genes for each cell identity
colLabels(sce) <- colData(sce)$CellType.Final
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

# Only keep genes that are relevant for each cell identity
mgs_fil <- lapply(names(mgs), function(i) {
    cat(i, "\n")
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.5, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)

# Down-sample sce data (to speed up computation)
down_sample_flag <- TRUE
# split cell indices by identity
if(down_sample_flag){
idx <- split(seq(ncol(sce)), sce$CellType.Final)
# downsample to at most 20 per identity & subset
n_cells <- 1000
cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
        n_cells <- n
    sample(i, n_cells)
})
sce <- sce[, unlist(cs_keep)]}

# Run SPOTLight
# Learn reference model
mod_ls <- trainNMF(
    x = sce,
    y = spe,
    groups = as.character(sce$CellType.Final),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

saveRDS(mod_ls, file.path(out_path, "NMF_model.rds"))

res <- runDeconvolution(
    x = spe,
    mod = mod_ls[["mod"]],
    ref = mod_ls[["topic"]])

saveRDS(res, file.path(out_path, "deconvolution_results.rds"))