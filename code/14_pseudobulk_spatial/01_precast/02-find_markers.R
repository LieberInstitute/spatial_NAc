library(here)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(getopt)
library(edgeR)
library(scran)
library(scater)
library(dplyr)
library(PCAtools)
library(gridExtra)
library(ggforce)
library(pheatmap)
library(pals)

spec <- matrix(
    c("agg_level", "a", 1, "character", "Aggregating at the donor or capture area level?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

cluster_col = paste0('precast_k', k)
ensembl_col = 'gene_id'
symbol_col = 'gene_name'
if(opt$agg_level == "sample_id_original"){
    spe_pseudo_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_capture_area', sprintf('spe_pseudo_%s.rds', cluster_col))
    modeling_rdata_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_capture_area', sprintf('model_results_%s.Rdata', cluster_col))
    modeling_genes_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_capture_area', sprintf('model_results_%s_FDR5perc.csv', cluster_col))
    plot_dir = here('plots', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_capture_area', cluster_col)
}else{
    if(opt$agg_level == "sample_id"){
        spe_pseudo_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_donor', sprintf('spe_pseudo_%s.rds', cluster_col))
        modeling_rdata_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_donor', sprintf('model_results_%s.Rdata', cluster_col))
        modeling_genes_path = here(
    'processed-data', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_donor', sprintf('model_results_%s_FDR5perc.csv', cluster_col))
        plot_dir = here('plots', '14_pseudobulk_spatial', '01_precast', 'pseudobulk_donor', cluster_col)
    }else{
        stop("Incorrect aggregation level")
    }
}

spe_pseudo = readRDS(spe_pseudo_path)
colData(spe_pseudo) <- colData(spe_pseudo)[ , unique(colnames(colData(spe_pseudo)))]
spe_pseudo$sample_id_original <- as.character(spe_pseudo$sample_id_original)
spe_pseudo$slide_num <- as.character(spe_pseudo$slide_num)
spe_pseudo[[cluster_col]] <- as.character(spe_pseudo[[cluster_col]])
spe_pseudo$sample_id_original <- make.names(spe_pseudo$sample_id_original)
spe_pseudo$slide_num <- make.names(spe_pseudo$slide_num)
spe_pseudo$sample_id_original <- factor(spe_pseudo$sample_id_original)
spe_pseudo$slide_num <- factor(spe_pseudo$slide_num)
spe_pseudo[[cluster_col]] <- make.names(spe_pseudo[[cluster_col]])
spe_pseudo[[cluster_col]] <- factor(spe_pseudo[[cluster_col]], levels = paste0("X", c(1:k)))

pdf(file = file.path(plot_dir, "histogram_boxplot_domain.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200, xlab ="Number of spots in each pseudobulk sample", ylab = "Frequency", main = "Distribution of # spots in each pseudobulk sample")
boxplot(ncells ~ colData(spe_pseudo)[[cluster_col]], data = colData(spe_pseudo), xlab = "Cluster", ylab = "Number of spots")
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id_original)
rowData(spe_pseudo)$high_expr_group_domain <- filterByExpr(spe_pseudo, group = colData(spe_pseudo)[[cluster_col]])

summary(rowData(spe_pseudo)$high_expr_group_sample_id)

summary(rowData(spe_pseudo)$high_expr_group_domain)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_domain))
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_sample_id, ]
dim(spe_pseudo)

x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x

spe_pseudo <- scuttle::addPerCellQC(
    spe_pseudo,
    subsets = list(Mito = which(grepl("^MT-", rowData(spe_pseudo)$gene_name)))
)
spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=3))
spe_pseudo<-spe_pseudo[,spe_pseudo$det_out==F]
dim(spe_pseudo)
rm(x)

#run PCA
set.seed(12141)
#runPCA(spe_pseudo)
pca <- prcomp(t(assays(spe_pseudo)$logcounts))

message(Sys.time(), " % of variance explained by the PCs:")
metadata(spe_pseudo)
metadata(spe_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca))
metadata(spe_pseudo)

pca_pseudo <- pca$x
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(spe_pseudo) <- list(PCA = pca_pseudo)


pdf(file = file.path(plot_dir, "pseudobulk_PCA_visium_final.pdf"), width = 14, height = 14)
plotPCA(spe_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sample_id_original", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = cluster_col, ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

# Extract the PCA data from spe_pseudo into a dataframe
PCAData <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
PCAData$domain <- colData(spe_pseudo)[[cluster_col]]

pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
    geom_point(aes(color = domain),size=2) + # Assuming you have a "domain" column in your PCAData
    labs(x = "PC1", y = "PC2") +
    scale_color_manual(values=unname(glasbey())[1:k], breaks = levels(PCAData$domain)) 
    theme_bw()+theme(panel.grid.major = element_blank(),  # Remove major grid lines
                     panel.grid.minor = element_blank())  # Remove minor grid lines

pdf(file= file.path(plot_dir,'pca_plot_final_check.pdf'),h=6,w=8)
pca_plot + geom_mark_ellipse(aes(color = domain, label = domain),
                             expand = unit(0.5,"mm"),
                             label.buffer = unit(-5, 'mm'),
                             show.legend = FALSE,
                             con.type = 'none')
dev.off()

outlier_clusters <- names(table(spe_pseudo[[cluster_col]]))[table(spe_pseudo[[cluster_col]]) < 10] 
spe_pseudo <- spe_pseudo[ ,!as.character(spe_pseudo[[cluster_col]]) %in% outlier_clusters]
spe_pseudo[[cluster_col]] <- factor(as.character(spe_pseudo[[cluster_col]]))
# Run DE model
# Rescale the feature age
spe_pseudo$age_scaled <- scales::rescale(spe_pseudo$Age,to=c(0,1))
if(!is.factor(spe_pseudo$Sex)){
    spe_pseudo$Sex <- factor(spe_pseudo$Sex, levels = c("F", "M"))
}
if(!is.factor(spe_pseudo$slide_num)){
    spe_pseudo$slide_num <- factor(spe_pseudo$slide_num, levels = unique(spe_pseudo$slide_num))
}
if(!is.factor(spe_pseudo$donor)){
    spe_pseudo$donor <- factor(spe_pseudo$donor, levels = unique(spe_pseudo$donor))
}

if(opt$agg_level == "sample_id_original"){
    covars <- c('Sex', 'slide_num')
}

if(opt$agg_level == "sample_id"){
    covars <- c('Sex', 'age_scaled')
}


mod<-registration_model(spe_pseudo, covars = covars, var_registration = cluster_col)
cors<-registration_block_cor(spe_pseudo, registration_model = mod,  var_sample_id = opt$agg_level)

results_enrichment <-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = cluster_col,
    var_sample_id = opt$agg_level,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

if(k >= 3){
    results_anova <-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = cluster_col,
    var_sample_id = opt$agg_level,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col, 
    suffix = "all"
)
}else{
    results_anova <- NULL
}

results_pairwise <-registration_stats_pairwise(
    spe_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = cluster_col,
    var_sample_id = opt$agg_level,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)

modeling_results <- list(
    "anova" = results_anova,
    "enrichment" = results_enrichment,
    "pairwise" = results_pairwise
)

save(modeling_results, file = modeling_rdata_path)

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
