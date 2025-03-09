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

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control")

# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]
# Read in the reference data set
sce_ext <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
sce_ext_neurons <- subset(sce_ext,idents = c("Drd1-MSN-1","Drd1-MSN-2","Drd2-MSN-1","Drd2-MSN-2",
                                           "Drd3-MSN","Grm8-MSN","GABAergic","Chat-Interneuron",
                                           "Pvalb-Interneuron","Sst-Interneuron","Glutamatergic"))
#Change the hyphen to a period
Idents(sce_ext_neurons) <- gsub(Idents(sce_ext_neurons),pattern = "-",replacement = ".")
Idents(sce_ext)   <- gsub(Idents(sce_ext),pattern = "-",replacement = ".")

# Processing the human registration results
reg_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
# Load snRNA-seq registration results
sn_registration <- readRDS(file.path(reg_dir, "sn_cellType_registration.rds"))
sn_enrichment <- sn_registration[["enrichment"]]
if(opt$subset_neurons){
  select_cell_types <- c("DRD1_MSN_A", "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D", "DRD2_MSN_A", 
"DRD2_MSN_B", "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F")
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
if(opt$subset_neurons){
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext_neurons), ]
}else{
  DEGs_df <- DEGs_df[DEGs_df$Rat_Gene_Symbol %in% rownames(sce_ext), ]
}

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


# Compute the t-stats for the rodent data
if(opt$subset_neurons){
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

