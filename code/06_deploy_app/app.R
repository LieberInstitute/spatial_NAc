library("spatialLIBD")
library("markdown")
library("here")

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds')

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
spe = readRDS(spe_path)
vars = colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
    spe,
    title = "Spatial NAc",
    spe_discrete_vars = c(
        vars[grep("^scran_", vars)],
        vars[grep("^10x_", vars)]
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    default_cluster = "path_groups",
    docs_path = "www"
)
