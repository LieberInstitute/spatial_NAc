####code for obtaining gene set enrichment for MCP factors
library(dplyr)
library(ggplot2)
library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(Seurat)
library(RColorBrewer)
library(spatialLIBD)
library(jaffelab)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(getopt)
library(edgeR)
library(scran)
library(scuttle)

set.seed(123)

# Read in the data
dat_dir <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(dat_dir)

plot_dir <- here::here("plots", "25_revisions","01_MCP_GSEA")
save_dir <- here::here("processed-data", "25_revisions", "01_MCP_GSEA")
dir.create(plot_dir, recursive = TRUE)
dir.create(save_dir, recursive = TRUE)

# Prepare expression data
expr_mat <- logcounts(spe)
geneData <- rowData(spe)
rownames(expr_mat) <- geneData$gene_id
# Only retain unique, protein coding autosomal genes
geneData <- geneData[!duplicated(geneData$gene_name), ]
geneData <- geneData[geneData$gene_type == "protein_coding", ]
geneData <- geneData[!grepl("^MT-", geneData$gene_name), ]
geneData <- geneData[!grepl("^RP[SL]", geneData$gene_name), ]

# Subset the expression matrix
expr_mat <- expr_mat[rownames(expr_mat) %in% geneData$gene_id, ]

# Now compute the correlations
pattern_cols <- grep("^meringue_cluster_", colnames(colData(spe)), value = TRUE)

compute_cor_matrix <- function(expr, patterns) {
  cor_mat <- matrix(NA, nrow = nrow(expr), ncol = length(patterns))
  rownames(cor_mat) <- rownames(expr)
  colnames(cor_mat) <- patterns

  for (pat in patterns) {
    print(pat)
    pat_values <- colData(spe)[[pat]]
    non_na_idx <- which(!is.na(pat_values))
    expr_sub <- expr[, non_na_idx, drop = FALSE]
    pat_sub <- pat_values[non_na_idx]
    expr_sub <- as.matrix(expr_sub)
    cors <- apply(expr_sub, 1, function(x) cor(x, pat_sub, method = "pearson"))
    cor_mat[, pat] <- cors}
  cor_mat
}

cor_mat <- compute_cor_matrix(expr_mat, pattern_cols)

cor_mat <- cor_mat[rowSums(is.na(cor_mat)) == 0, ]

# Prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h

################# GSEA
non0.mcp1 = rownames(cor_mat)[cor_mat[,"meringue_cluster_1"]>0]
non0.1.id = mapIds(org.Hs.eg.db, keys=non0.mcp1, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.1.id) = non0.mcp1
non0.1.id = non0.1.id[!is.na(non0.1.id)]
pathways.1 <- reactomePathways(non0.1.id)
pathways.1 <- x.h[names(pathways.1)]

mcp1.stats = cor_mat[names(non0.1.id),"meringue_cluster_1"]
names(mcp1.stats) = non0.1.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
mcp1.results = fgseaMultilevel(pathways.1, stats=mcp1.stats, scoreType="pos", minSize=15, maxSize=500)
mcp1.results$leadingEdge2 = sapply(mcp1.results$leadingEdge, paste, collapse="/")
write.csv(mcp1.results[,c(1:7,9)], file.path(save_dir ,"mcp1_reactome_results.csv"))

# Repeat for MCP2
non0.mcp2 = rownames(cor_mat)[cor_mat[,"meringue_cluster_2"]>0]
non0.2.id = mapIds(org.Hs.eg.db, keys=non0.mcp2, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.2.id) = non0.mcp2
non0.2.id = non0.2.id[!is.na(non0.2.id)]
pathways.2 <- reactomePathways(non0.2.id)
pathways.2 <- x.h[names(pathways.2)]

mcp2.stats = cor_mat[names(non0.2.id),"meringue_cluster_2"]
names(mcp2.stats) = non0.2.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
mcp2.results = fgseaMultilevel(pathways.2, stats=mcp2.stats, scoreType="pos", minSize=15, maxSize=500)
mcp2.results$leadingEdge2 = sapply(mcp2.results$leadingEdge, paste, collapse="/")
write.csv(mcp2.results[,c(1:7,9)], file.path(save_dir ,"mcp2_reactome_results.csv"))

# Repeat for MCP3
non0.mcp3 = rownames(cor_mat)[cor_mat[,"meringue_cluster_3"]>0]
non0.3.id = mapIds(org.Hs.eg.db, keys=non0.mcp3, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.3.id) = non0.mcp3
non0.3.id = non0.3.id[!is.na(non0.3.id)]
pathways.3 <- reactomePathways(non0.3.id)
pathways.3 <- x.h[names(pathways.3)]

mcp3.stats = cor_mat[names(non0.3.id),"meringue_cluster_3"]
names(mcp3.stats) = non0.3.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
mcp3.results = fgseaMultilevel(pathways.3, stats=mcp3.stats, scoreType="pos", minSize=15, maxSize=500)
mcp3.results$leadingEdge2 = sapply(mcp3.results$leadingEdge, paste, collapse="/")
write.csv(mcp3.results[,c(1:7,9)], file.path(save_dir ,"mcp3_reactome_results.csv"))

# Repeat for MCP4
non0.mcp4 = rownames(cor_mat)[cor_mat[,"meringue_cluster_4"]>0]
non0.4.id = mapIds(org.Hs.eg.db, keys=non0.mcp4, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
names(non0.4.id) = non0.mcp4
non0.4.id = non0.4.id[!is.na(non0.4.id)]
pathways.4 <- reactomePathways(non0.4.id)
pathways.4 <- x.h[names(pathways.4)]

mcp4.stats = cor_mat[names(non0.4.id),"meringue_cluster_4"]
names(mcp4.stats) = non0.4.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
mcp4.results = fgseaMultilevel(pathways.4, stats=mcp4.stats, scoreType="pos", minSize=15, maxSize=500)
mcp4.results$leadingEdge2 = sapply(mcp4.results$leadingEdge, paste, collapse="/")
write.csv(mcp4.results[,c(1:7,9)], file.path(save_dir ,"mcp4_reactome_results.csv"))

