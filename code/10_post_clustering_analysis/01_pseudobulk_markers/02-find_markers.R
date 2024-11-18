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
    c("cluster_col", "c", 1, "character", "Specify the name of the cluster column", 
      "pseudo_path", "p", 1, "character", "Specify the path of the pseudobulked expression",  
      "agg_level", "a", 1, "character", "Aggregating at the donor or capture area level?"),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

print(opt$cluster_col)
print(opt$pseudo_path)
print(opt$agg_level)
#opt <- list()
#opt$cluster_col <- 'BayesSpace_harmony_k12'
#opt$pseudo_path <- '02_BayesSpace/pseudobulk_capture_area'
#opt$agg_level <- "sample_id_original"
ensembl_col = 'gene_id'
symbol_col = 'gene_name'

if(grepl("BayesSpace", opt$cluster_col)){
    clust_col_n <- as.numeric(gsub("BayesSpace_harmony_k", "", opt$cluster_col))
    clust_col_n_nice <- sprintf("%02d", clust_col_n)
    opt$cluster_col <- paste0("BayesSpace_harmony_k", clust_col_n_nice)
}


spe_pseudo_path = here(
    'processed-data', '10_post_clustering_analysis', '01_pseudobulk_markers', 
     opt$pseudo_path, sprintf('spe_pseudo_%s.rds', opt$cluster_col))
modeling_rdata_path = here(
    'processed-data', '10_post_clustering_analysis', '01_pseudobulk_markers', 
     opt$pseudo_path, sprintf('model_results_%s.Rdata', opt$cluster_col))
modeling_genes_path = here(
    'processed-data', '10_post_clustering_analysis','01_pseudobulk_markers',
     opt$pseudo_path, sprintf('model_results_%s_FDR5perc.csv', opt$cluster_col))
plot_dir = here('plots', '10_post_clustering_analysis', '01_pseudobulk_markers', 
                 opt$pseudo_path, opt$cluster_col)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
spe_pseudo = readRDS(spe_pseudo_path)
colData(spe_pseudo) <- colData(spe_pseudo)[ , unique(colnames(colData(spe_pseudo)))]
if(opt$agg_level == "sample_id_original"){
    spe_pseudo$sample_id_original <- as.character(spe_pseudo$sample_id_original)
    spe_pseudo$sample_id_original <- make.names(spe_pseudo$sample_id_original)
    spe_pseudo$sample_id_original <- factor(spe_pseudo$sample_id_original)
}
spe_pseudo$slide_num <- as.character(spe_pseudo$slide_num)
spe_pseudo$slide_num <- make.names(spe_pseudo$slide_num)
spe_pseudo$slide_num <- factor(spe_pseudo$slide_num)


if(grepl("BayesSpace", opt$pseudo_path)){
    nClusters <- as.numeric(gsub("BayesSpace_harmony_k", "", opt$cluster_col))
    opt$cluster_col <- "clusters"
}

spe_pseudo[[opt$cluster_col]] <- as.character(spe_pseudo[[opt$cluster_col]])
if(!grepl("BayesSpace", opt$pseudo_path)){
    if(grepl("final_clusters", opt$pseudo_path)){
    spe_pseudo[[opt$cluster_col]] <- make.names(spe_pseudo[[opt$cluster_col]])
    spe_pseudo[[opt$cluster_col]] <- factor(spe_pseudo[[opt$cluster_col]])
}else{
    spe_pseudo[[opt$cluster_col]] <- make.names(spe_pseudo[[opt$cluster_col]])
    spe_pseudo[[opt$cluster_col]] <- factor(spe_pseudo[[opt$cluster_col]], levels = paste0("X", c(1:as.numeric(gsub("precast_k", "", opt$cluster_col)))))
}
}else{
    spe_pseudo[[opt$cluster_col]] <- make.names(spe_pseudo[[opt$cluster_col]])
    spe_pseudo[[opt$cluster_col]] <- factor(spe_pseudo[[opt$cluster_col]], levels = paste0("X", c(1:nClusters)))
}


pdf(file = file.path(plot_dir, "histogram_boxplot_domain.pdf"), width = 10, height = 6)
hist(spe_pseudo$ncells, breaks = 200, xlab ="Number of spots in each pseudobulk sample", ylab = "Frequency", main = "Distribution of # spots in each pseudobulk sample")
par(mar = c(10, 5, 2, 2))
boxplot(ncells ~ colData(spe_pseudo)[[opt$cluster_col]], data = colData(spe_pseudo), xlab = "", ylab = "Number of spots", las =2)
dev.off()

# Remove pseudobulked samples with too few spots
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
dim(spe_pseudo)

#find a good expression cutoff using edgeR::filterByExpr
if(opt$agg_level == "sample_id"){
    rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
}
if(opt$agg_level == "sample_id_original"){
    rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id_original)
}
rowData(spe_pseudo)$high_expr_group_domain <- filterByExpr(spe_pseudo, group = colData(spe_pseudo)[[opt$cluster_col]])

summary(rowData(spe_pseudo)$high_expr_group_sample_id)

summary(rowData(spe_pseudo)$high_expr_group_domain)

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_domain))
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_domain, ]
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

# Remove samples with too few detected genes
spe_pseudo <- scuttle::addPerCellQC(
    spe_pseudo,
    subsets = list(Mito = which(grepl("^MT-", rowData(spe_pseudo)$gene_name)))
)
#if(opt$agg_level == "sample_id_original"){
#    spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=2, batch = colData(spe_pseudo)[[opt$cluster_col]]))
#}
#if(opt$agg_level == "sample_id"){
#    spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=3, batch = colData(spe_pseudo)[[opt$cluster_col]]))
#}
spe_pseudo$det_out <- spe_pseudo$detected < 2000
#spe_pseudo$det_out<-as.logical(isOutlier(spe_pseudo$detected,type='lower',nmads=3, batch = colData(spe_pseudo)[[opt$cluster_col]]))
spe_pseudo<-spe_pseudo[,!spe_pseudo$det_out]
dim(spe_pseudo)
rm(x)

# Remove any clusters with fewer than 10 samples from the analysis
if(opt$agg_level == "sample_id_original"){
    outlier_clusters <- names(table(spe_pseudo[[opt$cluster_col]]))[table(spe_pseudo[[opt$cluster_col]]) < 10] 
}
if(opt$agg_level == "sample_id"){
    outlier_clusters <- names(table(spe_pseudo[[opt$cluster_col]]))[table(spe_pseudo[[opt$cluster_col]]) < 4] 
}
spe_pseudo <- spe_pseudo[ ,!as.character(spe_pseudo[[opt$cluster_col]]) %in% outlier_clusters]
spe_pseudo[[opt$cluster_col]] <- factor(as.character(spe_pseudo[[opt$cluster_col]]))

# Print final spe dimension
print("Final SPE dimensions")
dim(spe_pseudo)

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


pdf(file = file.path(plot_dir, "pseudobulk_PCA_visium_final.pdf"), width = 8, height = 8)
plotPCA(spe_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
plotPCA(spe_pseudo, colour_by = "ncells", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
if(opt$agg_level == "sample_id_original"){
    plotPCA(spe_pseudo, colour_by = "sample_id_original", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
}
if(opt$agg_level == "sample_id"){
    plotPCA(spe_pseudo, colour_by = "sample_id", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
}
plotPCA(spe_pseudo, colour_by = opt$cluster_col, ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(spe_pseudo)$PCA_var_explained)
dev.off()

# Plot the relationship between the PCs and known meta-data variables
if(opt$agg_level == "sample_id_original"){
    pdf(file = file.path(plot_dir, "explanatory_variables_PCA.pdf"), width = 8, height = 6)
plotExplanatoryPCs(spe_pseudo, variables = c("Age", opt$cluster_col, "sum", 
"detected", "slide_num", "subsets_Mito_percent", "sample_id"), npcs_to_plot = 10)
dev.off()
}

if(opt$agg_level == "sample_id"){
    pdf(file = file.path(plot_dir, "explanatory_variables_PCA.pdf"), width = 8, height = 6)
plotExplanatoryPCs(spe_pseudo, variables = c("Age","Sex", opt$cluster_col, "sum", "detected", "slide_num", "subsets_Mito_percent"), npcs_to_plot = 10)
dev.off()
}


# Extract the PCA data from spe_pseudo into a dataframe
PCAData <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
PCAData$domain <- colData(spe_pseudo)[[opt$cluster_col]]

pca_plot <- ggplot(PCAData, aes(PC01, PC02)) +
    geom_point(aes(color = domain),size=2) + # Assuming you have a "domain" column in your PCAData
    labs(x = "PC1", y = "PC2") +
    scale_color_manual(values=unname(glasbey())[1:length(unique(spe_pseudo[[opt$cluster_col]]))], breaks = levels(PCAData$domain)) 
    theme_bw()+theme(panel.grid.major = element_blank(),  # Remove major grid lines
                     panel.grid.minor = element_blank())  # Remove minor grid lines

pdf(file= file.path(plot_dir,'pca_plot_final_check.pdf'),h=6,w=8)
pca_plot + geom_mark_ellipse(aes(color = domain),
                             expand = unit(0.5,"mm"),
                             show.legend = FALSE,
                             con.type = 'none')
dev.off()

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


mod<-registration_model(spe_pseudo, covars = covars, var_registration = opt$cluster_col)
cors<-registration_block_cor(spe_pseudo, registration_model = mod,  var_sample_id = opt$agg_level)

results_enrichment <-registration_stats_enrichment(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = opt$cluster_col,
    var_sample_id = opt$agg_level,
    gene_ensembl = ensembl_col,
    gene_name = symbol_col
)
k <- length(unique(spe_pseudo[[opt$cluster_col]]))
if(k >= 3){
    results_anova <-registration_stats_anova(
    spe_pseudo,
    block_cor=cors,
    covars = covars,
    var_registration = opt$cluster_col,
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
    var_registration = opt$cluster_col,
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


spe_pseudo$spatialLIBD <- spe_pseudo[[opt$cluster_col]]
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
