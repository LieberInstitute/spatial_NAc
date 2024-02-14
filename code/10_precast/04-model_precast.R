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

registration_mod <- registration_model(spe_pseudo, covars = covars)

block_cor <- registration_block_cor(
    spe_pseudo, registration_model = registration_mod
)

results_enrichment <- registration_stats_enrichment(
    spe_pseudo,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

results_pairwise <- registration_stats_pairwise(
    spe_pseudo,
    registration_model = registration_mod,
    block_cor = block_cor,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

#   ANOVA is only defined and meaningful for 3 or more groups
if (opt$k >= 3) {
    results_anova <- registration_stats_anova(
        spe_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = ensembl_col,
        gene_name = symbol_col
        suffix = "all"
    )
} else {
    results_anova <- NULL
}

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)

save(modeling_results, modeling_path)

session_info()
