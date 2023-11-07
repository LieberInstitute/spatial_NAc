library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(scran)
library(BiocParallel)
library(spatialNAcUtils)
library(scran)
library(scater)

sample_info_path = here('raw-data', 'sample_key_spatial_NAc.csv')
sample_info_path2 = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
transformed_dir = here('processed-data', '04_VisiumStitcher')
processed_dir = here("processed-data", "05_harmony_BayesSpace")
raw_in_path = here('processed-data', '05_harmony_BayesSpace', 'spe_raw.rds')
filtered_out_path = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
plot_dir = here('plots', '05_harmony_BayesSpace')
num_cores = 4
num_red_dims = 20

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

#   Compute outlier spots by library size
spe$scran_low_lib_size <-
    factor(
        isOutlier(
            spe$sum_umi,
            type = "lower",
            log = TRUE,
            batch = spe$sample_id_original
        ),
        levels = c("TRUE", "FALSE")
    )

plot_list = list()
for (donor in unique(spe$sample_id)) {
    plot_list[[donor]] = spot_plot(
        spe,
        sample_id = donor,
        title = donor,
        var_name = "scran_low_lib_size",
        include_legend = TRUE,
        is_discrete = TRUE
    )
}
pdf(file.path(plot_dir, "sample_aware_low_lib_size.pdf"))
print(plot_list)
dev.off()

#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    (colSums(assays(spe)$counts) > 0)
]

message("Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id_original,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(
    spe,
    clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))

message("Running logNormCounts()")
spe <- logNormCounts(spe)

################################################################################
#   Compute PCA
################################################################################

message("Running modelGeneVar()")
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

message("Running getTopHVGs()")
top.hvgs <- getTopHVGs(dec, prop = 0.1)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))

save(
    top.hvgs, top.hvgs.fdr5, top.hvgs.fdr1,
    file = file.path(processed_dir, "top.hvgs.Rdata")
)

message("Running runPCA()")
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs, ncomponents = num_red_dims)
Sys.time()

#   Plot variance explained
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
pdf(file.path(dir_plots, "pca_elbow.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
dev.off()

message("Adding PCA to 'reducedDims(spe)'")
pca = prcomp(t(assays(spe)$counts), rank = num_red_dims)
reducedDims(spe)$PCA = pca$x
metadata(spe)$PCA_var_explained = jaffelab::getPcaVars(pca)[1:num_red_dims]
message(
    sprintf("Percent of variance explained for the top %s PCs:", num_red_dims)
)
print(metadata(spe)$PCA_var_explained)

message("Adding GLM-PCA to 'reducedDims(spe)'")
reducedDims(spe)$GLMPCA = glmpca(
    assays(spe)$counts, num_red_dims, fam = "poi", minibatch = "stochastic"
)$factors

message("Saving filtered spe")
saveRDS(spe, filtered_out_path)

session_info()
