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

code_dir = here("code", "05_harmony_BayesSpace")
source(paste0(code_dir, "/localOutlier_2.R"))

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
#   Read in the data and add additional QC metrics, followed by filtering data
################################################################################
spe = readRDS(raw_in_path)
cat("Initial number of spots:", dim(spe)[2], "\n")

# Add the number of estimated cells per spot using Vistoseg
sample_info_path = here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv')
sample_info = read.csv(sample_info_path) |> as_tibble()

counts_list = list()
for (i in 1:nrow(sample_info)) {
    counts_list[[i]] = read.csv(
                file.path(sample_info$spaceranger_dir[i],
                'tissue_spot_counts.csv')
            ) |>
            as_tibble() |>
            mutate(sample_id = sample_info$sample_id[i])
}

counts = do.call(rbind, counts_list) |>
    mutate(key = paste(barcode, sample_id, sep = '_')) |>
    filter(key %in% spe$key)
counts <- data.frame(counts)
counts <- counts[match(spe$key, counts$key), ]
spe$Nmask_dark_blue <- counts$Nmask_dark_blue
spe$Pmask_dark_blue <- counts$Pmask_dark_blue
spe$CNmask_dark_blue <- counts$CNmask_dark_blue

# Subset individuals when testing code to save memory and time
subset_individuals <- FALSE
if(subset_individuals){
    spe <- spe[ ,spe$sample_id %in% c("Br3942", "Br8325", "Br6423")]
    spe$sample_id_original <- as.factor(as.character(spe$sample_id_original))
    spe$sample_id <- as.factor(as.character(spe$sample_id))
    spe$slide_num <- as.factor(as.character(spe$slide_num))
    spe$donor <- as.factor(as.character(spe$donor))
    spe$array_num <- as.factor(as.character(spe$array_num))
}

# Preliminary QC
# Spot QC
# Select spots that have non-zero gene expression and found in tissue
spe <- spe[, (colSums(assays(spe)$counts) > 0) & spe$in_tissue]
cat("Number of spots after preliminary QC:", dim(spe)[2], "\n")
## Metrics QC
metrics_qc <- function(spe) {

    qc_df <- data.frame(
        log2sum = log2(spe$sum_umi),
        log2detected = log2(spe$sum_gene),
        subsets_Mito_percent = spe$expr_chrM_ratio*100,
        sample_id = spe$sample_id_original
    )

    qcfilter <- DataFrame(
        low_lib_size = isOutlier(qc_df$log2sum, type = "lower", log = TRUE, batch = qc_df$sample_id),
        low_n_features = isOutlier(qc_df$log2detected, type = "lower", log = TRUE, batch = qc_df$sample_id),
        high_subsets_Mito_percent = isOutlier(qc_df$subsets_Mito_percent, type = "higher", batch = qc_df$sample_id)
    )
    qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent


    spe$scran_low_lib_size_low_mito <- factor(qcfilter$low_lib_size & qc_df$subsets_Mito_percent < 0.5, levels = c("TRUE", "FALSE"))


    spe$scran_discard <-
        factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
    spe$scran_low_lib_size <-
        factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
    spe$scran_low_n_features <-
        factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
    spe$scran_high_subsets_Mito_percent <-
        factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

    ## Find edge spots
    spots <- data.frame(
        row = spe$array_row,
        col = spe$array_col,
        sample_id = spe$sample_id_original
    )

    edge_spots_row <- group_by(spots, sample_id, row) %>% summarize(min_col = min(col), max_col = max(col))
    edge_spots_col <- group_by(spots, sample_id, col) %>% summarize(min_row = min(row), max_row = max(row))

    spots <- left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
    spots$edge_spots <- with(spots, row == min_row | row == max_row | col == min_col | col == max_col)

    spots$row_distance <- with(spots, pmin(abs(row - min_row), abs(row - max_row)))
    spots$col_distance <- with(spots, pmin(abs(col - min_col), abs(col - max_col)))
    ## spots$edge_distance <- with(spots, sqrt(row_distance^2 + col_distance^2))
    ## The above is from:
    ## sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
    ## but it was wrong, here's a case the the smallest distance is on the column:
    ## sqrt(0^2 + col_distance^2) = col_distance
    spots$edge_distance <- with(spots, pmin(row_distance, col_distance))


    spe$edge_spots <- factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
    spe$edge_distance <- spots$edge_distance


    spe$scran_low_lib_size_edge <- factor(qcfilter$low_lib_size & spots$edge_distance < 1, levels = c("TRUE", "FALSE"))

    return(spe)
}

spe <- metrics_qc(spe)
# Remove spots that are both edge spots and have a low library size
spe <- spe[ , !(spe$scran_low_lib_size == "TRUE" & spe$edge_spots == "TRUE")]

# Visualize QC metrics prior to any further filtering
pdf(width = 14, height = 5, paste0(plot_dir, "/QC_metric_prior_filtering.pdf"))
ggplot(data = as.data.frame(colData(spe)),
       aes(x = sum_umi)) +
  geom_histogram(aes(y = after_stat(density)), 
                 colour = "black", 
                 fill = "grey", bins = 50) +
  geom_density(alpha = 0.5,
               adjust = 1.0,
               fill = "#A0CBE8",
               colour = "#4E79A7") +
  geom_vline(xintercept = max(spe$sum_umi[spe$scran_low_lib_size == "TRUE"]), col="red", linetype = "dashed")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Library size") + 
  ylab("Density") + 
  theme_classic()  + theme(plot.margin = margin(2, 2, 2, 2, "cm"))

p <- ggplot(as.data.frame(colData(spe)), aes(x=Nmask_dark_blue, y=sum_umi)) +
  geom_point(size=0.5) + 
  geom_smooth(method = "loess", se=FALSE) +
  geom_hline(yintercept = 200, colour='red') + xlab("Nuclei count") + ylab("Sum UMI") +
  theme_minimal() 
ggMarginal(p, type='histogram', margins = 'both')

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

dev.off()

filterSpotSweeper <- FALSE
if(filterSpotSweeper){
# Check which spots are found to be outliers based on SpotSweeper
spe$array_row <- spe$array_row_original
spe$array_col <- spe$array_col_original
spe$sample_id_original <- as.character(spe$sample_id_original)
# Find SpotSweeper outliers
feats <- c("sum_umi", "sum_gene", "expr_chrM_ratio")
spe <- localOutliers_2(spe, features = feats, n_neighbors=18, 
                    data_output=TRUE,
                    method="multivariate", samples = "sample_id_original", n_cores = num_cores)
# Check if spotSweeper outliers are enriched for white matter

spe$array_row <- spe$array_row_transformed
spe$array_col <- spe$array_col_transformed
}

# Gene QC
rownames(spe) <- rowData(spe)$gene_name
exprLogic <- counts(spe) > 0
nCells_by_slide <- lapply(levels(spe$slide_num), function(iSlide){
    rowSums(exprLogic[ ,spe$slide_num == iSlide])
})
nCells_by_slide <- do.call(cbind, nCells_by_slide)
colnames(nCells_by_slide) <- levels(spe$slide_num)
select.genes <-  rownames(nCells_by_slide)[rowSums(nCells_by_slide == 0) == 0]
# Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]

# Compute the relative expression of each gene per cell
# Use sparse matrix operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
C <- counts(spe)
C@x <- C@x / rep.int(colSums(C), diff(C@p)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = .1, las = 1, xlab = "% total count per cell", col = scales::hue_pal()(20)[20:1], horizontal = TRUE)

# Exclude mitochondrial genes from downstream analysis
spe <- spe[!grepl("^MT-", rownames(spe)), ]
# Filter Ribossomal gene (optional if that is a problem on your data)
spe <- spe[!grepl("^RP[SL]", rownames(spe)), ]
spe <- spe[!rownames(spe) == "MALAT1", ]
# Filter Hemoglobin gene  (optional if that is a problem on your data)
#spe <- spe[!grepl("^HB[^(P|E|S)]", rownames(spe)), ]

################################################################################
#   Compute log-normalized counts
################################################################################
C <- counts(spe)
gene_attr <- data.frame(mean = rowMeans(C), detection_rate = rowMeans(C > 
    0), var = apply(C, 1, var))
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(C)
cell_attr <- data.frame(n_umi = colSums(C), n_gene = colSums(C > 
    0))
rownames(cell_attr) <- colnames(C)

ggplot(gene_attr, aes(log_mean, log_var)) + geom_point(alpha = 0.3, shape = 16) + 
    geom_density_2d(size = 0.3) + geom_abline(intercept = 0, slope = 1, color = "red")

#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
message(Sys.time(), " - Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$slide_num,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message(Sys.time(), " - Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(
    spe,
    clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores), 
    positive = FALSE
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
# We want to obtain genes that capture biologically interesting trends
# as opposed to genes that vary between batches. Since the donor and the 
# slide number are not perfectly correlated, we'll take the more conservative approach
# to look for highly variable genes in each capture area and summarize across
# the capture areas. 
dec <- modelGeneVar(
    spe,
    block = spe$donor,
    BPPARAM = MulticoreParam(num_cores), 
    density.weights = FALSE
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


################################################################################
#   Compute GLM-PCA
################################################################################

message(Sys.time(), " - Running devianceFeatureSelection()")
spe <- devianceFeatureSelection(
    spe, assay = "counts", fam = "binomial", sorted = FALSE, 
    batch = as.factor(spe$sample_id_original))

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

message(Sys.time(), " - Running nullResiduals()")
spe <- nullResiduals(
    # default params
    spe, assay = "counts", fam = "binomial", type = "deviance"
)

message(Sys.time(), " - Running GLM-PCA")
spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.2000, ncomponents = num_red_dims,
    name = "GLMPCA_approx",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.5000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_5000",
    BSPARAM = BiocSingular::IrlbaParam()
)

spe <- runPCA(
    spe,
    exprs_values = "binomial_deviance_residuals",
    subset_row = hdgs.hb.10000, ncomponents = num_red_dims,
    name = "GLMPCA_approx_10000",
    BSPARAM = BiocSingular::IrlbaParam()
)

################################################################################
#   Obtain preliminary clusters based on default GLM-PCA and PCA settings
################################################################################
spe$leiden20_PCA <- clusterCells(spe, 
                           use.dimred = "PCA", 
                           BLUSPARAM = SNNGraphParam(k = 20, 
                                                     cluster.fun = "leiden"), 
                                                     cluster.args = list(resolution=2.0))

spe$leiden20_GLMPCA <- clusterCells(spe, 
                           use.dimred = "GLMPCA_approx", 
                           BLUSPARAM = SNNGraphParam(k = 20, 
                                                     cluster.fun = "leiden", 
                                                     cluster.args = list(resolution=2.0)))


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
