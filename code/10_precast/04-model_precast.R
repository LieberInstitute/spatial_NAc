library(here)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(getopt)
library(spatialLIBD)

spec <- matrix(
    c("k", "k", 1, "integer", "Number of PRECAST clusters"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

cluster_col = paste0('precast_k', opt$k)
covars = c("Age", "Sex")
ensembl_col = 'gene_id'
symbol_col = 'gene_name'

spe_pseudo_path = here(
    'processed-data', '10_precast', sprintf('spe_pseudo_%s.rds', cluster_col)
)
modeling_path = here(
    'processed-data', '10_precast',
    sprintf('modeling_results_%s.Rdata', cluster_col)
)

spe_pseudo = readRDS(spe_pseudo_path)

#   Modeling code

save(modeling_results, modeling_path)

session_info()
