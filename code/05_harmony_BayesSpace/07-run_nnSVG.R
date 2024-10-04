library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)
library(nnSVG)
library(scran)
library(scuttle)
library(getopt)

spec <- matrix(
    c(
        "use_precast", "p", 1, "logical", "Use k=2 PRECAST results as a covariate to nnSVG?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
message(Sys.time(), " | Loading SpatialExperiment")
spe <- loadHDF5SummarizedExperiment(spe_dir)
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

if (opt$use_precast) {
    out_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_precast_out",
        paste0(sample_id, ".csv")
    )
    precast_path <- here("processed-data", "10_precast", "00_pre_clustering", "PRECAST_k2.csv")
} else {
    out_path <- here(
        "processed-data", "05_harmony_BayesSpace", "07-run_nnSVG", "nnSVG_out",
        paste0(sample_id, ".csv")
    )
}

set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   Subset to this sample, drop overlapping spots, and bring into memory
#-------------------------------------------------------------------------------

message(Sys.time(), " | Subsetting to this sample and bringing into memory")
spe <- spe[, (spe$sample_id == sample_id) & !spe$exclude_overlapping]

#   Bring fully into memory to speed up computations later
assays(spe) <- list(counts = assays(spe)$counts)
spe <- realize(spe)

#-------------------------------------------------------------------------------
#   Import PRECAST k=2 results if applicable, and only take spots classified by
#   PRECAST
#-------------------------------------------------------------------------------

if (opt$use_precast) {
    message("Will use PRECAST k=2 clusters as a covariate to nnSVG.")

    precast_results <- read.csv(precast_path) |>
        as_tibble()

    #   Drop any spots that don't have a PRECAST cluster defined
    defined_spots <- spe$key %in% precast_results$key
    message(
        sprintf(
            "%s | Dropping %s spots, which did not have a PRECAST cluster defined.",
            Sys.time(),
            length(which(!defined_spots))
        )
    )
    spe <- spe[, defined_spots]
}

#-------------------------------------------------------------------------------
#   Filter lowly expressed and mitochondrial genes, and take spots with at least
#   some nonzero counts
#-------------------------------------------------------------------------------

message(Sys.time(), " | Filtering genes and spots")
spe <- filter_genes(
    spe,
    filter_genes_ncounts = 3,
    filter_genes_pcspots = 0.5,
    filter_mito = TRUE
)
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]
message("Dimensions of spe after filtering:")
print(dim(spe))

#-------------------------------------------------------------------------------
#   Recompute logcounts (library-size normalization as recommended in
#   https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html)
#-------------------------------------------------------------------------------

message(Sys.time(), " | Re-computing logcounts")
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

#-------------------------------------------------------------------------------
#   Run nnSVG and export results
#-------------------------------------------------------------------------------

message(Sys.time(), " | Running nnSVG")
if (opt$use_precast) {
    precast_cluster <- colData(spe) |>
        as_tibble() |>
        left_join(precast_results, by = "key") |>
        pull(cluster)
    precast_cluster_ <- as.factor(precast_cluster)
    capture_area_ <- as.factor(as.character(spe$sample_id_original))
    spe <- nnSVG(spe, X = model.matrix(~precast_cluster_ + capture_area_))
} else {
    capture_area_ <- as.factor(as.character(spe$sample_id_original))
    spe <- nnSVG(spe, X = model.matrix(~capture_area_))
}

message(Sys.time(), " | Exporting results")
write.csv(rowData(spe), out_path, row.names = FALSE, quote = FALSE)

session_info()

