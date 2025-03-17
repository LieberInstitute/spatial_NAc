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
        "data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)

opt <- getopt(spec)
#opt <- list()
#opt$data  <- "human_NAc"

print(opt$data)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "cross_validation", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "cross_validation", opt$data)

sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

print("The dimensions of the sce object")
print(dim(sce))

if(opt$data == "human_NAc"){
    test_ranks <- c(5,10,20, 30, 40, 50, 75, 100,125)
}else{
    test_ranks <- c(5, 10, 20, 30, 40, 50, 60, 70, 80)
}


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
    test_density = 0.2)

dir.create(res_dir)
saveRDS(cvnmf, file = file.path(res_dir,paste0("nmf_cv_results.rds")))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(plot_dir, paste0("nmf_cv_results.pdf")))
plot(cvnmf)
dev.off()

