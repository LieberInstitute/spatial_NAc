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

spec <- matrix(
    c(
        "subset_neurons", "n", 1, "logical", "Perform correlation only for neuronal cell types?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
#opt <- list()
#opt$subset_neurons <- FALSE

print(opt$subset_neurons)
dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "non_human_primate")
plot_dir <- here("plots", "17_cross_species_comparison", "non_human_primate")

# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]
# Read in the reference data set
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")

if(opt$subset_neurons){
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
}else{
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}
if(opt$subset_neurons){
  pdf(file.path(plot_dir, "QC_metrics_NHP_neurons.pdf"), width = 5, height = 4)
  VlnPlot(sce_ext, features = "nFeature_RNA", group.by = "MSN_type", pt.size = 0) + ggtitle("Number of genes")+ theme(legend.position = "none") + xlab("")
  dev.off()
}else{
  pdf(file.path(plot_dir, "QC_metrics_NHP.pdf"), width = 7, height = 4)
  VlnPlot(sce_ext, features = "nCount.RNA", group.by = "cell_type", pt.size = 0) + ggtitle("UMI count") + theme(legend.position = "none") + xlab("")
  VlnPlot(sce_ext, features = "nFeatures_RNA", group.by = "cell_type", pt.size = 0) + ggtitle("Number of genes") + theme(legend.position = "none") + xlab("")
  dev.off()
}

#Change the hyphen to a period
if(opt$subset_neurons){
   Idents(sce_ext) <- sce_ext$MSN_type
   Idents(sce_ext)   <- gsub(Idents(sce_ext),pattern = "-",replacement = ".")
}else{
  Idents(sce_ext) <- sce_ext$cell_type
  Idents(sce_ext)   <- gsub(Idents(sce_ext),pattern = "-",replacement = ".")
}

# Processing the human registration results
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
# Load snRNA-seq registration results
sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]
if(opt$subset_neurons){
  select_cell_types <- c("DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A", 
"DRD2_MSN_B")
}else{
  select_cell_types <- c("DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A", 
"DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F","Astrocyte_A", "Astrocyte_B",
 "Microglia", "OPC", "Oligo", "Ependymal", "Endothelial")
}

sce <- sce[ ,sce$CellType.Final %in% select_cell_types]
DEG_Lists <- list()

for(i in select_cell_types){
  df <- sn_enrichment[ ,grep(i, colnames(sn_enrichment))]
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

geneData <- rowData(sce)
DEGs_df$Gene_symbol <- NA
for(i in 1:nrow(DEGs_df)){
 DEGs_df$Gene_symbol[i] <- geneData$gene_name[geneData$gene_id == DEGs_df$Gene_ID[i]]
}
# Add in the NHP symbols
Monkey_Orthos <- orthologs(genes = DEGs_df$Gene_symbol,species = "Macaca mulatta")
colnames(Monkey_Orthos)[5] <- "Monkey_Symbol"

#Now merge the DEGs_df and the monkey symbol. 
DEGs_df <-merge(x = DEGs_df,
                            y = Monkey_Orthos,
                            by.x = "Gene_symbol",
                            by.y = "human_symbol")

DEGs_df <- na.omit(DEGs_df)
# Remove any non one-to-one apping
select_human_genes <- names(table(DEGs_df$Gene_symbol)[table(DEGs_df$Gene_symbol) == 1])
select_nhp_genes <- names(table(DEGs_df$Monkey_Symbol)[table(DEGs_df$Monkey_Symbol)==1])
DEGs_df <- DEGs_df[DEGs_df$Gene_symbol %in% select_human_genes, ]
DEGs_df <- DEGs_df[DEGs_df$Monkey_Symbol %in% select_nhp_genes, ]
# Only select those genes that we have expression in the rat data for
DEGs_df <- DEGs_df[DEGs_df$Monkey_Symbol %in% rownames(sce_ext[["RNA"]]$data), ]


rownames(sce) <- rowData(sce)$gene_name

t_stat_mat_human <- matrix(nrow = length(unique(DEGs_df$Gene_symbol)),ncol = length(select_cell_types))
#Change the row names and col names
rownames(t_stat_mat_human) <- DEGs_df$Gene_symbol
colnames(t_stat_mat_human) <- select_cell_types
for(i in rownames(t_stat_mat_human)){
  print(i)
  for(l in colnames(t_stat_mat_human)){
    #First get all other celltypes
    All_Others <- select_cell_types[which(select_cell_types != l)]
    #Pull gene expression values for cells in cluster of interest
    c1_vals <- logcounts(sce)[i,sce$CellType.Final == l]
    #Pull gene expression values of all other cells
    c2_vals <- logcounts(sce)[i,sce$CellType.Final %in% All_Others]
    #Get the mean for cell type of interest
    x_c1  <- mean(c1_vals)
    #Get the mean for all other cell types
    x_c2 <- mean(c2_vals)
    #Get standard deviation of c1 vals
    sdev_c1 <- sd(x = c1_vals)
    #Get standard deviation of c1c cals
    sdec_c2 <- sd(x = c2_vals)
    #Now calculate t
    t_stat_mat_human[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(sum(sce$CellType.Final == l))
  }
}
if(opt$subset_neurons){
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human_neurons.rds"))
}else{
  saveRDS(t_stat_mat_human, file.path(res_dir, "t_stat_mat_human.rds"))
}


# Compute the t-stats for the NHP data
t_stat_mat_nhp <- matrix(nrow = length(unique(DEGs_df$Monkey_Symbol)),ncol = length(levels(sce_ext)))
#Change the row names and col names
rownames(t_stat_mat_nhp) <- DEGs_df$Monkey_Symbol
colnames(t_stat_mat_nhp) <- levels(sce_ext)
for(i in rownames(t_stat_mat_nhp)){
  print(i)
  for(l in colnames(t_stat_mat_nhp)){
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
    t_stat_mat_nhp[i,l] <- (x_c1-x_c2)/sqrt((sdev_c1^2+sdec_c2^2)/2)*sqrt(length(WhichCells(sce_ext,idents = l)))
  }
}
if(opt$subset_neurons){
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp_neurons.rds"))
}else{
  saveRDS(t_stat_mat_nhp, file.path(res_dir, "t_stat_mat_nhp.rds"))
}
