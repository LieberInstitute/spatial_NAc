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
qc_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', '02-compute_QC_metrics', 'spe_with_QC_metrics_hdf5'
)
filtered_ordinary_path = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe','spe_filtered.rds'
)
filtered_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe','spe_filtered_hdf5'
)
plot_dir = here('plots', '05_harmony_BayesSpace', '03-filter_normalize_spe')
num_red_dims = 50

num_cores = Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(0)

################################################################################
#   Read in the data and add additional QC metrics, followed by filtering data
################################################################################
spe = loadHDF5SummarizedExperiment(qc_hdf5_dir)
cat("Initial number of spots:", dim(spe)[2], "\n")

unique_samples <- unique(spe$sample_id)
pdf(width = 8, height = 8, paste0(plot_dir, "/local_outliers.pdf"))
spot_plot(spe, "Br2720", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br2743", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br3942", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6423", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6432", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6471", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6522", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8325", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8492", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8667", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)
dev.off()

# Visualize QC metrics prior to any further filtering
pdf(width = 12, height = 5, paste0(plot_dir, "/QC_metrics_prior_filtering.pdf"))
ggplot(data = as.data.frame(colData(spe)),
       aes(x = sample_id_original, y = sum_umi, fill = sample_id_original)) + 
  geom_violin() +
  xlab("Capture Area") + 
  ylab("Library size") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Capture Area"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = donor, y = sum_umi, fill = donor)) + 
  geom_violin() +
  xlab("Donor") + 
  ylab("Library size") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Donor"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = slide_num, y = sum_umi, fill = slide_num)) + 
  geom_violin() +
  xlab("Slide Number") + 
  ylab("Library size") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Slide Number"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = sample_id_original, y = sum_gene, fill = sample_id_original)) + 
  geom_violin() +
  xlab("Capture Area") + 
  ylab("# Detected genes") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Capture Area"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = donor, y = sum_gene, fill = donor)) + 
  geom_violin() +
  xlab("Donor") + 
  ylab("# Detected genes") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Donor"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = slide_num, y = sum_gene, fill = slide_num)) + 
  geom_violin() +
  xlab("Slide Number") + 
  ylab("# Detected genes") + 
  theme_classic()  + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2, "cm")) +
  guides(fill=guide_legend(title="Slide Number"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = sum_umi)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey", bins = 50) +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  geom_vline(xintercept = 200, col="red", linetype = "dashed")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()  + theme(plot.margin = margin(2, 2, 2, 2, "cm"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = sum_gene)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey", bins = 50) +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Genes expressed in each spot") + 
  ylab("Density") + 
  theme_classic()  + theme(plot.margin = margin(2, 2, 2, 2, "cm"))

ggplot(data = as.data.frame(colData(spe)),
       aes(x = expr_chrM_ratio)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey", bins = 50) +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Ratio of mitochondrial expression") + 
  ylab("Density") + 
  theme_classic() + theme(plot.margin = margin(2, 2, 2, 2, "cm"))

p <- ggplot(as.data.frame(colData(spe)), aes(x=Nmask_dark_blue, y=sum_umi)) +
  geom_point(size=0.5) + 
  geom_smooth(method = "loess", se=FALSE) +
  geom_hline(yintercept = 200, colour='red') + xlab("Nuclei count") + ylab("Sum UMI") +
  theme_bw() 

p2 <- ggplot(as.data.frame(colData(spe)), aes(x=Nmask_dark_blue, y=sum_gene)) +
  geom_point(size=0.5) + 
  geom_smooth(method = "loess", se=FALSE) +
  geom_hline(yintercept = 200, colour='red') + xlab("Nuclei count") + ylab("Number of detected genes") +
  theme_bw() 

p3 <- ggplot(as.data.frame(colData(spe)), aes(x=Nmask_dark_blue, y=expr_chrM_ratio)) +
  geom_point(size=0.5) + 
  geom_smooth(method = "loess", se=FALSE) + xlab("Nuclei count") + ylab("Ratio of mitochondrial gene expression") +
  theme_bw() 

plot_grid(ggMarginal(p, type='histogram', margins = 'both'), nrow = 1, ncol = 1)
plot_grid(ggMarginal(p2, type='histogram', margins = 'both'), nrow = 1, ncol = 1)
plot_grid(ggMarginal(p3, type='histogram', margins = 'both'), nrow = 1, ncol = 1)

dev.off()

spe$low_umi <- spe$sum_umi < 250
spe$low_gene <- spe$sum_gene < 250
spe$low_gene_edge_spot <- spe$low_umi & spe$edge_distance < 6

pdf(width = 8, height = 8, paste0(plot_dir, "/spots_low_ngenes.pdf"))
spot_plot(spe, "Br2720", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br2743", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br3942", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6423", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6432", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6471", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6522", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8325", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8492", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8667", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE)
dev.off()

pdf(width = 8, height = 8, paste0(plot_dir, "/spots_low_umi.pdf"))
spot_plot(spe, "Br2720", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br2743", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br3942", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6423", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6432", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6471", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6522", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8325", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8492", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8667", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE)
dev.off()

# Visualize edge spots that are low library size and would be exckuded
pdf(width = 8, height = 8, paste0(plot_dir, "/edge_spots_low_lib_size.pdf"))
spot_plot(spe, "Br2720", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br2743", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br3942", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6423", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6432", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6471", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br6522", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8325", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8492", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
spot_plot(spe, "Br8667", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE)
dev.off()

# Perform optional additional spot filtering
spe <- spe[ ,spe$scran_low_lib_size_edge == "FALSE"]
# Check if SpotSweeper outliers are enriched for white matter
lost <- calculateAverage(counts(spe)[,spe$local_outliers == "TRUE"])
kept <- calculateAverage(counts(spe)[,!spe$local_outliers == "TRUE"])
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

is_WM <- rowData(spe)$gene_name %in% c("MBP", "MOBP", "SLC47A2", "SERPINA5")
pdf(width = 5, height = 5, paste0(plot_dir, "/diff_expression_spotSweeper_QC.pdf"))
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[is_WM], logFC[is_WM], col="dodgerblue", pch=16)
dev.off()
# Since spotSweeper doesn't seem to be selectively picking White matter, we will exclude these spots
filterSpotSweeper <- TRUE
if(filterSpotSweeper){
spe <- spe[ ,spe$local_outliers == "FALSE"]
}

# Gene QC
exprLogic <- counts(spe) > 0
spe$sample_id <- as.factor(spe$sample_id)
nSpots_by_donor <- lapply(levels(spe$sample_id), function(iSample){
    rowSums(exprLogic[ ,spe$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]

# Compute the relative expression of each gene per cell
# Use sparse matrix operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
C <- as(counts(spe, withDimnames = TRUE), "dgCMatrix")
C@x <- C@x / rep.int(colSums(C), diff(C@p)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

pdf(width = 8, height = 8, paste0(plot_dir, "/most_expressed_genes.pdf"))
par(mar = c(3, 8, 3, 3))
boxplot(as.matrix(t(C[most_expressed, ])), cex = .1, las = 1, xlab = "% total count per cell", col = scales::hue_pal()(20)[20:1], horizontal = TRUE)
dev.off()
# Exclude mitochondrial genes from downstream analysis
spe <- spe[!grepl("^MT-", rowData(spe)$gene_name), ]
# Filter Ribossomal gene (optional if that is a problem on your data)
spe <- spe[!grepl("^RP[SL]", rowData(spe)$gene_name), ]
spe <- spe[!rowData(spe)$gene_name == "MALAT1", ]
# Filter Hemoglobin gene  (optional if that is a problem on your data)
#spe <- spe[!grepl("^HB[^(P|E|S)]", rownames(spe)), ]


################################################################################
#   Compute log-normalized counts
################################################################################
C <- counts(spe, withDimnames = TRUE)
mean_expr <- rowMeans(C)
det_rate <- rowMeans(C > 0)
var_expr <- rowVars(C)
gene_attr <- data.frame(mean = mean_expr, detection_rate = det_rate, var = var_expr)
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- make.names(rownames(C), unique = TRUE)
cell_attr <- data.frame(n_umi = colSums(C), n_gene = colSums(C > 
    0))
rownames(cell_attr) <- colnames(C)

pdf(width = 8, height = 8, paste0(plot_dir, "/mean_var_counts.pdf"))
ggplot(gene_attr, aes(log_mean, log_var)) + geom_point(alpha = 0.3, shape = 16) + 
    geom_density_2d(size = 0.3) + geom_abline(intercept = 0, slope = 1, color = "red")
dev.off()

#message(Sys.time(), " - Running quickCluster()")

#Sys.time()
#spe$scran_quick_cluster <- quickCluster(
#    spe,
#    BPPARAM = MulticoreParam(num_cores),
#    block = spe$slide_num,
#    block.BPPARAM = MulticoreParam(num_cores)
#)
Sys.time()

message(Sys.time(), " - Running computeLibraryFactors()")
Sys.time()
#spe <- computeSumFactors(
#    spe,
#    clusters = spe$scran_quick_cluster,
#    BPPARAM = MulticoreParam(num_cores), 
#    positive = FALSE
#)
spe <- computeLibraryFactors(spe)
Sys.time()

#print("Quick cluster table:")
#table(spe$scran_quick_cluster)

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
#   Visualize spots by slide number and donor
################################################################################
donor <- c()
slide <- c()
nSpots <- c()
spe$slide_num <- as.character(spe$slide_num)
for(d in levels(spe$donor)){
    for(s in unique(spe$slide_num[spe$donor == d])){
        donor <- c(donor, d)
        slide <- c(slide, s)
        nSpots <- c(nSpots, sum(spe$donor == d & spe$slide_num == s))
    }
}
spot_summary <- data.frame("donor" = donor, "slide" = slide, "nSpots" = nSpots)
spe$slide_num <- as.factor(spe$slide_num)

pdf(file.path(plot_dir, "spots_by_donor_after_filtering.pdf"), useDingbats = FALSE)
ggplot(spot_summary, aes(fill=slide, y=nSpots, x=donor)) + 
    geom_bar(position="stack", stat="identity") + theme_classic() + scale_fill_jco() + xlab("Donor") + ylab("Number of spots")
dev.off()

################################################################################
#   Compute PCA
################################################################################

message(Sys.time(), " - Running modelGeneVar()")
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
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
pdf(file.path(plot_dir, "pca_elbow.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p1"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p1.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p1", ylab = "Variance explained (%)")
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p2"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p2.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p2", ylab = "Variance explained (%)")
dev.off()

percent.var <- attr(reducedDim(spe, "PCA_p5"), "percentVar")
pdf(file.path(plot_dir, "pca_elbow_p5.pdf"), useDingbats = FALSE)
plot(percent.var, xlab = "PC_p5", ylab = "Variance explained (%)")
dev.off()

pdf(file.path(plot_dir, "pca_explanatory_vars.pdf"), useDingbats = FALSE)
plotExplanatoryPCs(spe,variables = c("sum_umi", "sum_gene", "donor", "slide_num", 
                                    "Nmask_dark_blue"), npcs_to_plot = 10, dimred = "PCA_p2") 
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

#message(Sys.time(), " - Running nullResiduals()")
#spe <- nullResiduals(
    # default params
#    spe, assay = "counts", fam = "binomial", type = "deviance",
#    batch = as.factor(spe$slide_num)
#)

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000, ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)
message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.5000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_5000",
    BSPARAM = BiocSingular::IrlbaParam()
)
message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.10000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_10000",
    BSPARAM = BiocSingular::IrlbaParam()
)

################################################################################
#   Save the processed SPE object
################################################################################

message(Sys.time(), " - Saving HDF5-backed filtered spe")
spe = saveHDF5SummarizedExperiment(
    spe, dir = filtered_hdf5_dir, replace = TRUE
)

message(Sys.time(), " - Saving ordinary filtered spe")
spe = realize(spe)
saveRDS(spe, filtered_ordinary_path)

################################################################################
#   Examine reduced dimensions
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

pdf(file.path(plot_dir, "GLMPCA_approx.pdf"), width = 6, height = 6)
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "donor")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "slide_num")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "sum_umi")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "MOBP")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "low_umi")
plotReducedDim(spe, dimred = "GLMPCA_approx", ncomponents = 2, colour_by = "low_gene")
dev.off()

session_info()


