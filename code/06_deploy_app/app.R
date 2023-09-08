library("spatialLIBD")
library("markdown")

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
spe <- readRDS("spe_shiny.rds")
vars <- colnames(colData(spe))
spe$exclude_overlapping <- as.factor(spe$exclude_overlapping)

## Deploy the website
spatialLIBD::run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "Spatial NAc",
    spe_discrete_vars = c(
        vars[grep("^scran_", vars)],
        vars[grep("^10x_", vars)],
        "exclude_overlapping",
        "ManualAnnotation"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    default_cluster = "10x_graphclust",
    auto_crop_default = FALSE,
    docs_path = "www"
)
