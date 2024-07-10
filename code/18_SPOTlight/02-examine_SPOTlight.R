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

# Read in the results from SPOTlight
mod_ls <- readRDS(file.path(out_path, "NMF_model.rds"))
res <- readRDS(file.path(out_path, "deconvolution_results.rds"))


