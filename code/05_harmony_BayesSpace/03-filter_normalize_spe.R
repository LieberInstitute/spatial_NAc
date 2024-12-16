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
spot_plot(spe, "Br2720", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE) +
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br2720)")
spot_plot(spe, "Br2743", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br2743)")
spot_plot(spe, "Br3942", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br3942)")
spot_plot(spe, "Br6423", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br6423)")
spot_plot(spe, "Br6432", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br6432)")
spot_plot(spe, "Br6471", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br6471)")
spot_plot(spe, "Br6522", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br6522)")
spot_plot(spe, "Br8325", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br8325)")
spot_plot(spe, "Br8492", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br8492)")
spot_plot(spe, "Br8667", var_name = "local_outliers", is_discrete = TRUE, spatial = TRUE)+
            guides(fill = guide_legend(override.aes = list(size = 5))) + ggtitle("Local outliers (Br8667)")
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
spot_plot(spe, "Br2720", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2720") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br2743", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2743") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br3942", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br3942") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6423", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6423") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6432", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6432") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6471", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6471") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6522", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6552") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8325", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8325") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8492", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8492") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8667", var_name = "low_gene", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8667") + guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

pdf(width = 8, height = 8, paste0(plot_dir, "/spots_low_umi.pdf"))
spot_plot(spe, "Br2720", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2720") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br2743", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2743") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br3942", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br3942") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6423", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6423") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6432", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6432") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6471", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6471") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6522", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6522") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8325", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8325") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8492", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8492") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8667", var_name = "low_umi", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8667") + guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

# Visualize edge spots that are low library size and would be exckuded
pdf(width = 8, height = 8, paste0(plot_dir, "/edge_spots_low_lib_size.pdf"))
spot_plot(spe, "Br2720", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2720") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br2743", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br2743") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br3942", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br3942") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6423", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6423") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6432", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6432") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6471", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6471") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br6522", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br6522") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8325", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8325") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8492", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8492") + guides(fill = guide_legend(override.aes = list(size = 5)))
spot_plot(spe, "Br8667", var_name = "low_gene_edge_spot", is_discrete = TRUE, spatial = TRUE) + ggtitle("Br8667") + guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

# Perform optional additional spot filtering
spe <- spe[ ,spe$low_gene_edge_spot == "FALSE"]
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
rownames(C) <- rowData(spe)$gene_name
C@x <- C@x / rep.int(colSums(C), diff(C@p)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

pdf(width = 8, height = 8, paste0(plot_dir, "/most_expressed_genes.pdf"))
par(mar = c(3, 8, 3, 3))
boxplot(as.matrix(t(C[most_expressed, ])), cex = .1, las = 1, xlab = "% total count per cell", col = scales::hue_pal()(20)[20:1], horizontal = TRUE)
dev.off()
# Exclude mitochondrial genes from downstream analysis
#spe <- spe[!grepl("^MT-", rowData(spe)$gene_name), ]
# Filter Ribossomal gene (optional if that is a problem on your data)
#spe <- spe[!grepl("^RP[SL]", rowData(spe)$gene_name), ]
#spe <- spe[!rowData(spe)$gene_name == "MALAT1", ]
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

session_info()


