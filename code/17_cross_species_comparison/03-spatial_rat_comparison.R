library(SpatialExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(Seurat)
library(RColorBrewer)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(getopt)
library(babelgene)
library(pheatmap)
library(HDF5Array)
library(ComplexHeatmap)

spe_dir <- here("processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5")
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control_spatial")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control_spatial")
subset_neurons <- FALSE

print(subset_neurons)
# Read in human 10x visium
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Read in external data
dat_dir <- here("processed-data", "12_snRNA")
sce_ext <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
sce_ext_neurons <- subset(sce_ext,idents = c("Drd1-MSN-1","Drd1-MSN-2","Drd2-MSN-1","Drd2-MSN-2",
                                           "Drd3-MSN","Grm8-MSN","GABAergic","Chat-Interneuron",
                                           "Pvalb-Interneuron","Sst-Interneuron","Glutamatergic"))
#Change the hyphen to a period
Idents(sce_ext_neurons) <- gsub(Idents(sce_ext_neurons),pattern = "-",replacement = ".")
Idents(sce_ext)   <- gsub(Idents(sce_ext),pattern = "-",replacement = ".")

# Processing the human registration results
reg_dir <- here("processed-data", "10_post_clustering_analysis","01_pseudobulk_markers", "01_precast", "pseudobulk_capture_area", "final_clusters")
# Load snRNA-seq registration results
load(file.path(reg_dir, "model_results_precast_clusters.Rdata"))
sp_enrichment <- modeling_results$enrichment
rm(modeling_results)

if(subset_neurons){
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", 
"Inhibitory")
}else{
  select_domains <- c("MSN.1", "MSN.2", "MSN.3", "D1.islands", "Excitatory", 
"Inhibitory", "WM", "Endothelial.Ependymal")
}

# Add precast results to the spe object
clusters_file <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["precast_clusters"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_file), by = 'key') |>
    pull(cluster) |>
    as.factor()
# Remove spots with no PRECAST output
spe <- spe[ ,!is.na(spe[["precast_clusters"]])]

# Replace any spatial characters with a .
spe[["precast_clusters"]] <- gsub(" ", ".", spe[["precast_clusters"]])
spe[["precast_clusters"]] <- gsub("/", ".", spe[["precast_clusters"]])

spe <- spe[ ,spe$precast_clusters %in% select_domains]
DEG_Lists <- list()

for(i in select_domains){
  df <- sp_enrichment[ ,grep(i, colnames(sp_enrichment))]
  colnames(df) <- gsub(paste0("_", i), "", colnames(df))
  df <- df[df$fdr < 0.05, ]
  df <- df[df$logFC > 0, ]
  df <- df[order(df$logFC, decreasing = TRUE), ]
  df <- df[1:min(250,dim(df)[1]) , ]
  df$CellType <- i
  df$Gene_ID <- rownames(df)
  DEG_Lists[[i]] <- df
}
#Create a dataframe
DEGs_df <- do.call(what = rbind, DEG_Lists)
rownames(DEGs_df) <- c(1:dim(DEGs_df)[1])
geneData <- rowData(spe)
DEGs_df$Gene_symbol <- NA
for(i in 1:nrow(DEGs_df)){
 DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]]
}
# Add in the mouse symbols
# Add human symbol corresponding to the rat gene
rat_human_ortholog <- read.delim("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/raw-data/HOM_AllOrganism.rpt")
hom_rat <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "rat", ]
hom_hs <- rat_human_ortholog[rat_human_ortholog$Common.Organism.Name == "human", ]

DEGs_df$Rat_Gene_Symbol <- NA
for(i in 1:nrow(DEGs_df)){
  print(i)
  #Only run the code if a homolog is found. 
  if(length(hom_rat[hom_rat$DB.Class.Key %in% hom_hs[hom_hs$Symbol %in% DEGs_df[i,"Gene_symbol"],"DB.Class.Key"],"Symbol"]) == 1){
    DEGs_df[i,"Rat_Gene_Symbol"] <- hom_rat[hom_rat$DB.Class.Key %in% hom_hs[hom_hs$Symbol %in% DEGs_df[i,"Gene_symbol"],"DB.Class.Key"],"Symbol"]
  }else{
    next
  }
}
DEGs_df <- na.omit(DEGs_df)

# Remove any non one-to-one apping
select_human_genes <- names(table(DEGs_df$Gene_symbol)[table(DEGs_df$Gene_symbol) == 1])
select_rat_genes <- names(table(DEGs_df$Rat_Gene_Symbol)[table(DEGs_df$Rat_Gene_Symbol)==1])
DEGs_df <- DEGs_df[DEGs_df$Gene_symbol %in% select_human_genes, ]
DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% select_rat_genes, ]
# Only select those genes that we have expression in the rat data for
if(subset_neurons){
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext_neurons), ]
}else{
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext), ]
}

rownames(spe) <- rowData(spe)$gene_name

t_stat_mat_human <- matrix(nrow = length(unique(DEGs_df$Gene_symbol)),ncol = length(select_domains))
#Change the row names and col names
rownames(t_stat_mat_human) <- DEGs_df$Gene_symbol
colnames(t_stat_mat_human) <- select_domains
for(i in rownames(t_stat_mat_human)){
  print(i)
  for(l in colnames(t_stat_mat_human)){
    #First get all other celltypes
    All_Others <- select_domains[which(select_domains != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- logcounts(spe)[i,spe$precast_clusters == l]
    #Pull gene expression values of all other cells
    c2_vals <- logcounts(spe)[i,spe$precast_clusters %in% All_Others]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_human[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(sum(spe$precast_clusters == l))
  }
}
if(subset_neurons){
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons.rds"))
}else{
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human.rds"))
}

# Compute the t-stats for the rodent data
if(subset_neurons){
  t_stat_mat_rat <- matrix(nrow = length(unique(DEGs_df$Rat_Gene_Symbol)),ncol = length(levels(sce_ext_neurons)))
  rownames(t_stat_mat_rat) <- DEGs_df$Rat_Gene_Symbol
  colnames(t_stat_mat_rat) <- levels(sce_ext_neurons)
  for(i in rownames(t_stat_mat_rat)){
  print(i)
  for(l in colnames(t_stat_mat_rat)){
    #First get all other celltypes
    All_Others <- levels(sce_ext_neurons)[which(levels(sce_ext_neurons) != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- GetAssayData(object = sce_ext_neurons,assay = "RNA",slot = "data")[i,WhichCells(sce_ext_neurons,idents = l)]
    #Pull gene expression values of all other cells
    c2_vals <- GetAssayData(object = sce_ext_neurons,assay = "RNA",slot = "data")[i,WhichCells(sce_ext_neurons,idents = All_Others)]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_rat[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(WhichCells(sce_ext_neurons,idents = l)))
  }
}
saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat_neurons.rds"))
}else{
  t_stat_mat_rat <- matrix(nrow = length(unique(DEGs_df$Rat_Gene_Symbol)),ncol = length(levels(sce_ext)))
  rownames(t_stat_mat_rat) <- DEGs_df$Rat_Gene_Symbol
  colnames(t_stat_mat_rat) <- levels(sce_ext)
  for(i in rownames(t_stat_mat_rat)){
  print(i)
  for(l in colnames(t_stat_mat_rat)){
    #First get all other celltypes
    All_Others <- levels(sce_ext)[which(levels(sce_ext) != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- GetAssayData(object = sce_ext,assay = "RNA",slot = "data")[i,WhichCells(sce_ext,idents = l)]
    #Pull gene expression values of all other cells
    c2_vals <- GetAssayData(object = sce_ext,assay = "RNA",slot = "data")[i,WhichCells(sce_ext,idents = All_Others)]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_rat[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(WhichCells(sce_ext,idents = l)))
  }
}
saveRDS(t_stat_mat_rat, file.path(res_dir, "t_stat_mat_rat.rds"))
}

# Plot 
t_stat_mat_rat <- t_stat_mat_rat[rowSums(is.na(t_stat_mat_rat)) == 0, ]
t_stat_mat_human <- t_stat_mat_human[rowSums(is.na(t_stat_mat_human)) == 0, ]

DEGs_df <- DEGs_df[DEGs_df$Gene_symbol %in% rownames(t_stat_mat_human), ]
DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(t_stat_mat_rat), ]

t_stat_mat_rat <- t_stat_mat_rat[rownames(t_stat_mat_rat) %in% DEGs_df$Rat_Gene_Symbol, ]
t_stat_mat_rat <- t_stat_mat_rat[match(DEGs_df$Rat_Gene_Symbol, rownames(t_stat_mat_rat)), ]

t_stat_mat_human <- t_stat_mat_human[rownames(t_stat_mat_human) %in% DEGs_df$Gene_symbol, ]
t_stat_mat_human <- t_stat_mat_human[match(DEGs_df$Gene_symbol, rownames(t_stat_mat_human)), ]

corr_mat <- matrix(NA, nrow = dim(t_stat_mat_rat)[2],ncol = dim(t_stat_mat_human)[2])
rownames(corr_mat) <- colnames(t_stat_mat_rat)
colnames(corr_mat) <- colnames(t_stat_mat_human)
for(i in rownames(corr_mat)){
  for(j in colnames(corr_mat)){
    corr_mat[i, j] <- cor(as.numeric(t_stat_mat_rat[ ,i]), as.numeric(t_stat_mat_human[ ,j]))
  }
}

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("\\.", " ", colnames(corr_mat))

if(subset_neurons){
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory")
  row_order <- c("Drd1 MSN 1", "Drd2 MSN 1", "Drd1 MSN 2",  "Drd2 MSN 2", "Drd3 MSN", "Grm8 MSN", "Glutamatergic", "GABAergic", "Pvalb Interneuron", "Chat Interneuron","Sst Interneuron")
}else{
  col_order <- c("MSN 1", "MSN 2", "MSN 3", "D1 islands", "Excitatory", "Inhibitory", "WM", "Endothelial Ependymal")
  row_order <- c("Drd1 MSN 1", "Drd2 MSN 1", "Drd1 MSN 2",  "Drd2 MSN 2", "Drd3 MSN", "Grm8 MSN",
   "Glutamatergic", "GABAergic", "Pvalb Interneuron", "Chat Interneuron","Sst Interneuron",
   "Astrocyte", "Microglia", "Polydendrocyte", "Olig 1", "Mural" )
}

corr_mat <- corr_mat[ ,match(col_order, colnames(corr_mat))]
corr_mat <- corr_mat[match(row_order, rownames(corr_mat)), ]

# Plot the heatmap of correlation
col_fun <- circlize::colorRamp2(c(-0.65,0,0.65),hcl_palette = "Purple-Green")
complex_plot_cor <- ComplexHeatmap::Heatmap(matrix = t(corr_mat),
                              name = "Correlation",
                              column_title = NULL,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_title = NULL,
                              rect_gp = gpar(col = "white", lwd = 2),
                              col = col_fun, 
                              border = TRUE,
                              heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(6, "cm"),
                              title_position = "topcenter", title_gp = gpar(fontsize = 14), border = "black", at = c(-0.65, -0.2, 0.2, 0.65)), 
                              row_names_side = "left")

if(subset_neurons){
  pdf(file.path(plot_dir, "correlation_heatmap_sn_neurons.pdf"), width = 5, height = 4)
  draw(complex_plot_cor, heatmap_legend_side="bottom")
  dev.off()
}else{
  pdf(file.path(plot_dir, "correlation_heatmap_sn.pdf"), width = 7, height = 5)
  draw(complex_plot_cor, heatmap_legend_side="bottom")
  dev.off()
}



