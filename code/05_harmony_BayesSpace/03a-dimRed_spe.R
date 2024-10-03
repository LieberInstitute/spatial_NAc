library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(scran)
library(scater)
library(scry)
library(HDF5Array)
library(bluster)
library(SpotSweeper)
library(pheatmap)
library(ggspavis)
library(ggplot2)
library(ggExtra)
library(edgeR)
library(cowplot)
library(ggsci)
library(dittoSeq)

processed_dir = here("processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe")
filtered_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe','spe_filtered_hdf5'
)
dimRed_ordinary_path = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe','spe_filtered_dimRed.rds'
)
dimRed_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe','spe_filtered_dimRed_hdf5'
)
num_red_dims = 50
plot_dir = here('plots', '05_harmony_BayesSpace', '03-filter_normalize_spe')

num_cores = Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(0)
################################################################################
#   Read in the data and add additional QC metrics, followed by filtering data
################################################################################
spe = loadHDF5SummarizedExperiment(filtered_hdf5_dir)
################################################################################
#   Visualize spots by slide number and donor
################################################################################
donor <- c()
capture_area <- c()
nSpots <- c()
spe$sample_id_original <- as.character(spe$sample_id_original)
for(d in levels(spe$donor)){
    for(s in unique(spe$sample_id_original[spe$donor == d])){
        donor <- c(donor, d)
        capture_area <- c(capture_area, s)
        nSpots <- c(nSpots, sum(spe$donor == d & spe$sample_id_original == s))
    }
}
spot_summary <- data.frame("donor" = donor, "capture_area" = capture_area, "nSpots" = nSpots)
spe$sample_id_original <- as.factor(spe$sample_id_original)

pdf(file.path(plot_dir, "spots_by_donor_after_filtering.pdf"), useDingbats = FALSE)
ggplot(spot_summary, aes(fill=capture_area, y=nSpots, x=donor)) + 
    geom_bar(position="dodge", stat="identity") + theme_classic() + xlab("Donor") + ylab("Number of spots")
dev.off()

################################################################################
#   Compute PCA
################################################################################

message(Sys.time(), " - Running modelGeneVar()")
# From
# http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
# We want to obtain genes that capture biologically interesting trends
# as opposed to genes that vary between batches. Since the donor and the 
# slide number are not perfectly correlated, we'll take the more conservative approach
# to look for highly variable genes in each capture area and summarize across
# the capture areas. 
dec <- modelGeneVar(
    spe,
    block = spe$sample_id,
    BPPARAM = MulticoreParam(num_cores), 
    density.weights = FALSE
)

rowData(spe) <- cbind(rowData(spe), dec)
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
top.hvgs.p1 <- getTopHVGs(dec, prop = 0.1)
top.hvgs.p2 <- getTopHVGs(dec, prop = 0.2)
top.hvgs.p5 <- getTopHVGs(dec, prop = 0.5)

top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
print(paste("Num HVGs at FDR = 0.05:", length(top.hvgs.fdr5)))

top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
print(paste("Num HVGs at FDR = 0.01:", length(top.hvgs.fdr1)))


save(
    top.hvgs.p1, top.hvgs.p2, top.hvgs.p5, top.hvgs.fdr5, top.hvgs.fdr1,
    file = file.path(processed_dir, "top.hvgs.Rdata")
)

message(Sys.time(), " - Running runPCA()")
Sys.time()
spe <- runPCA(spe, subset_row = top.hvgs.fdr5, ncomponents = num_red_dims, name = "PCA")
spe <- runPCA(spe, subset_row = top.hvgs.p1, ncomponents = num_red_dims, name = "PCA_p1")
spe <- runPCA(spe, subset_row = top.hvgs.p2, ncomponents = num_red_dims, name = "PCA_p2")
spe <- runPCA(spe, subset_row = top.hvgs.p5, ncomponents = num_red_dims, name = "PCA_p5")
Sys.time()

#   Plot variance explained
percent.var <- attr(reducedDim(spe, "PCA"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow.pdf"), useDingbats = FALSE, width = 5, height = 5)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)", pch = 20)
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p1"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p1.pdf"), useDingbats = FALSE, width = 5, height = 5)
plot(percent.var, xlab = "PC_p1", ylab = "Variance explained (%)", pch = 20)
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p2"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p2.pdf"), useDingbats = FALSE, width = 5, height = 5)
plot(percent.var, xlab = "PC_p2", ylab = "Variance explained (%)", pch = 20)
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p5"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p5.pdf"), useDingbats = FALSE, width = 5, height = 5)
plot(percent.var, xlab = "PC_p5", ylab = "Variance explained (%)", pch = 20)
dev.off()

pdf(file.path(plot_dir, "pca_explanatory_vars.pdf"), useDingbats = FALSE, width = 10, height = 5)
plotExplanatoryPCs(spe,variables = c("sample_id_original", "sum_umi", "sum_gene", "donor", "slide_num", 
                                    "Nmask_dark_blue"), npcs_to_plot = 10, dimred = "PCA", theme_size = 15)
dev.off()

################################################################################
#   Examine PCA
################################################################################
rbPal <- colorRampPalette(c('red','blue'))
rownames(spe) <- rowData(spe)$gene_name
cols <- rbPal(100)[as.numeric(cut(assays(spe)$logcounts["MBP", ], 100))]
pdf(file.path(plot_dir, "PC1_vs_UMI.pdf"), width = 5, height = 5)
plot(spe$sum_umi,reducedDim(spe, "PCA_p2")[ ,1], pch = 20,col = alpha(cols, 0.4), xlab = "Sum UMI", 
ylab = "PC1")
dev.off()

pdf(file.path(plot_dir, "PCA_p2.pdf"), width = 6, height = 6)
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "donor")
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "slide_num")
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "sum_umi")
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "MOBP")
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "low_umi")
plotReducedDim(spe, dimred = "PCA_p2", ncomponents = 2, colour_by = "low_gene")
dev.off()

################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")

spe <- devianceFeatureSelection(
    spe, assay = "counts", fam = "binomial", sorted = FALSE, 
    batch = as.factor(spe$sample_id))

pdf(file.path(plot_dir, "binomial_deviance.pdf"))
plot(sort(rowData(spe)$binomial_deviance, decreasing = T),
    type = "l", xlab = "ranked genes",
    ylab = "binomial deviance", main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()
hdgs.hb.2000 <- rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:2000]
hdgs.hb.5000 <- rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:5000]
hdgs.hb.10000 <- rownames(spe)[order(rowData(spe)$binomial_deviance, decreasing = TRUE)][1:10000]

save(
    hdgs.hb.2000, hdgs.hb.5000, hdgs.hb.10000,
    file = file.path(processed_dir, "hdgs.hb.Rdata")
)

residuals <- lapply(levels(spe$sample_id), function(i){
    cat(i, "\n")
    spe_subset <- subset(spe, ,sample_id == i)
    spe_subset <- nullResiduals(spe_subset, assay = "counts", fam = "binomial", type = "deviance")
    assay(spe_subset, "binomial_deviance_residuals")
})
residuals <- do.call(cbind, residuals)
residuals <- residuals[ ,match(colnames(spe), colnames(residuals))]
assay(spe, "binomial_deviance_residuals") <- residuals


message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000, ncomponents = num_red_dims,
    name = "GLMPCA_approx")
message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.5000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_5000")
message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.10000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_10000")

################################################################################
#   Save the processed SPE object
################################################################################

message(Sys.time(), " - Saving HDF5-backed filtered spe")
spe = saveHDF5SummarizedExperiment(
    spe, dir = dimRed_hdf5_dir, replace = TRUE
)

message(Sys.time(), " - Saving ordinary filtered spe")
spe = realize(spe)
saveRDS(spe, dimRed_ordinary_path)

################################################################################
#   Examine reduced dimensions
################################################################################
pdf(file.path(plot_dir, "GLMPCA_approx.pdf"), width = 6, height = 6)
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "donor")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "slide_num")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "sum_umi")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "MOBP")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "low_umi")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "low_gene")
dev.off()

session_info()


