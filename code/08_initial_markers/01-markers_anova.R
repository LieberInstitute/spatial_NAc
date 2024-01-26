library("here")
library("sessioninfo")
library("SpatialExperiment")
library("spatialLIBD")
library("tidyverse")

spe_in <- here("processed-data", "05_harmony_BayesSpace", "spe_harmony.rds")
cluster_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "clusters_BayesSpace"
)
k2_cluster_name <- "BayesSpace_harmony_k02"

#   Load SpatialExperiment and add k = 2 clusters
spe <- readRDS(spe_in)
cluster_import(
    spe,
    k2_cluster_name,
    cluster_dir = cluster_dir
)

covars <- c("BayesSpace", "age", "sex")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

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
results_anova <-
    registration_stats_anova(
        sce_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name,
        suffix = suffix
    )
