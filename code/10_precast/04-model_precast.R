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

registration_mod <-
    registration_model(sce_pseudo, covars = covars)

block_cor <-
    registration_block_cor(sce_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        sce_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
results_pairwise <-
    registration_stats_pairwise(
        sce_pseudo,
        registration_model = registration_mod,
        block_cor = block_cor,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
if (k >= 3) {
    results_anova <-
        registration_stats_anova(
            sce_pseudo,
            block_cor = block_cor,
            covars = covars,
            gene_ensembl = gene_ensembl,
            gene_name = gene_name,
            suffix = suffix
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
