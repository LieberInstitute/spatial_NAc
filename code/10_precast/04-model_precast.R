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
modeling_rdata_path = here(
    'processed-data', '10_precast',
    sprintf('model_results_%s.Rdata', cluster_col)
)
modeling_genes_path = here(
    'processed-data', '10_precast',
    sprintf('model_results_%s_FDR5perc.csv', cluster_col)
)

spe_pseudo = readRDS(spe_pseudo_path)

################################################################################
#   Use spatialLIBD registration functions to find interesting genes WRT
#   pairwise comparisons, enrichment, and ANOVA for gene expression
#   pseudobulked by PRECAST clusters
################################################################################

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

save(modeling_results, modeling_rdata_path)

################################################################################
#   Export CSV of significant genes
################################################################################

spe_pseudo$spatialLIBD <- spe_pseudo[[cluster_col]]
sig_genes <- sig_genes_extract_all(
    n = nrow(spe_pseudo),
    modeling_results = modeling_results,
    sce_layer = spe_pseudo
)

fix_csv <- function(df) {
    for (i in seq_len(ncol(df))) {
        if (any(grepl(",", df[, i]))) {
            message(paste(Sys.time(), "fixing column", colnames(df)[i]))
            df[, i] <- gsub(",", ";", df[, i])
        }
    }
    return(df)
}
z <- fix_csv(as.data.frame(subset(sig_genes, fdr < 0.05)))
write.csv(z[, !grepl("^in_rows", colnames(z))], file = modeling_genes_path)

session_info()
