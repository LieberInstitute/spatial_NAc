library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(scran)
library(scater)
library(scry)
library(HDF5Array)

processed_dir = here("processed-data", "05_harmony_BayesSpace")
raw_in_path = here('processed-data', '05_harmony_BayesSpace', 'spe_raw.rds')
filtered_ordinary_path = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
filtered_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
plot_dir = here('plots', '05_harmony_BayesSpace')
num_red_dims = 50

num_cores = Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(0)

################################################################################
#   Compute log-normalized counts
################################################################################

spe = readRDS(raw_in_path)

#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]

message(Sys.time(), " - Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id_original,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message(Sys.time(), " - Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(
    spe,
    clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

print("Quick cluster table:")
table(spe$scran_quick_cluster)

message(Sys.time(), " - Running checking sizeFactors()")
summary(sizeFactors(spe))

message(Sys.time(), " - Running logNormCounts()")
spe <- logNormCounts(spe)

#   Save a copy of the SPE with HDF5-backed assays, which will be important to
#   control memory consumption later
message(Sys.time(), " - Saving HDF5-backed object to control memory later")
spe = saveHDF5SummarizedExperiment(
    spe, dir = paste0(filtered_hdf5_dir, '_temp'), replace = TRUE
)
gc()

################################################################################
#   Compute PCA
################################################################################

message(Sys.time(), " - Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
dec <- modelGeneVar(
    spe,
    block = spe$sample_id_original,
    BPPARAM = MulticoreParam(num_cores)
)

pdf(file.path(plot_dir, "scran_modelGeneVar.pdf"), useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE
    )
}, dec$per.block, names(dec$per.block))
dev.off()

message(Sys.time(), " - Running getTopHVGs()")
top.hvgs <- getTopHVGs(dec, prop = 0.1)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))

save(
    top.hvgs, top.hvgs.fdr5, top.hvgs.fdr1,
    file = file.path(processed_dir, "top.hvgs.Rdata")
)

message(Sys.time(), " - Running runPCA()")
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = num_red_dims)
Sys.time()

#   Plot variance explained
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
dev.off()

################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(
    spe, assay = "counts", fam = "binomial", sorted = FALSE,
    batch = as.factor(spe$sample_id_original)
)

pdf(file.path(plot_dir, "binomial_deviance.pdf"))
plot(sort(rowData(spe)$binomial_deviance, decreasing = T),
    type = "l", xlab = "ranked genes",
    ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals(
    # default params
    spe, assay = "counts", fam = "binomial", type = "deviance"
)


hdgs.hb <- rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb, ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

message(Sys.time(), " - Saving HDF5-backed filtered spe")
spe = saveHDF5SummarizedExperiment(
    spe, dir = filtered_hdf5_dir, replace = TRUE
)

message(Sys.time(), " - Saving ordinary filtered spe")
spe = realize(spe)
saveRDS(spe, filtered_ordinary_path)

session_info()
