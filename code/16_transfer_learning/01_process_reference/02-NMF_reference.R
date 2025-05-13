library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(dplyr)
library(getopt)
library(Matrix)
library(HDF5Array)
set.seed(1234)

spec <- matrix(
    c("data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

#opt <- list()
#opt$data <- "human_NAc"
print(opt$data)

# Read data based on whether we are processing the human or rat NAc snRNA seq data
# Read data and create Seurat object
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)

sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

print("The dimensions of the sce object")
print(dim(sce))

dat <- sce[["RNA"]]$data

print("Running RCppML")
options(RcppML.threads=4)
if(opt$data == "human_NAc"){
  x <- RcppML::nmf(dat,
                 k=66,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

if(opt$data == "rat_case_control_acute"){
  x <- RcppML::nmf(dat,
                 k=43,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

if(opt$data == "rat_case_control_repeated"){
  x <- RcppML::nmf(dat,
                 k=40,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

if(opt$data == "rat_case_control_morphine_acute"){
   x <- RcppML::nmf(dat,
                 k=46,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

if(opt$data == "rat_case_control_morphine_repeated"){
   x <- RcppML::nmf(dat,
                 k=36,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 1135,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)
}

saveRDS(x, file = file.path(res_dir,paste0("nmf_results.rds")))
