library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(getopt)
library(ComplexHeatmap)
library(scater)
library(scran)
library(here)
library(scuttle)
library(Matrix)
library(SingleCellExperiment)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(RcppML)
library(singlet)
set.seed(1234)

spec <- matrix(
    c(
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs", 
        "data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)

opt <- getopt(spec)
#opt <- list()
#opt$gene_selection_strategy <- "all_genes"
#opt$data  <- "human_NAc"
#opt$replicate <- 1
print(opt$gene_selection_strategy)
print(opt$data)


# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "cross_validation", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "cross_validation", opt$data)


# Read in the spatial experiment data
raw_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5"
)
spe <- loadHDF5SummarizedExperiment(raw_in_path)
spe$low_umi <- spe$sum_umi < 250
spe$low_gene <- spe$sum_gene < 250
spe$low_gene_edge_spot <- spe$low_umi & spe$edge_distance < 6
spe <- spe[ ,spe$low_gene_edge_spot == "FALSE"]
spe <- spe[ ,spe$local_outliers == "FALSE"]
exprLogic <- counts(spe) > 0
spe$sample_id <- as.factor(spe$sample_id)
nSpots_by_donor <- lapply(levels(spe$sample_id), function(iSample){
    rowSums(exprLogic[ ,spe$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]


if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
  }else{
        stop("Invalid input data set")
  }
}

if(opt$gene_selection_strategy == "all_genes"){
  if(opt$data == "human_NAc"){
       sce <- sce[rowData(sce)$gene_id %in% rowData(spe)$gene_id, ]
  }
  if(opt$data == "rat_case_control"){
    refDir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
    orthologs_df <- readRDS(file.path(refDir, opt$data, "orthologs_df.rds"))
    sce <- sce[rownames(sce) %in% orthologs_df$rat_genes, ]
  }
}

test_ranks <- c(5,10,20, 30, 40, 50, 75, 100,125)
#test_ranks <- c(5, 10)
if(opt$data == "human_NAc"){
    cvnmf <- cross_validate_nmf(
    logcounts(sce),
    ranks=test_ranks,
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)
}

if(opt$data == "rat_case_control"){
    cvnmf <- cross_validate_nmf(
    sce[["RNA"]]$data,
    ranks=test_ranks,
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)
}

dir.create(res_dir)
saveRDS(cvnmf, file = file.path(res_dir,paste0("nmf_cv_results_", opt$gene_selection_strategy,".rds")))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plot_dir, paste0("nmf_cv_results_", opt$gene_selection_strategy, ".pdf")))
plot(cvnmf)
dev.off()

