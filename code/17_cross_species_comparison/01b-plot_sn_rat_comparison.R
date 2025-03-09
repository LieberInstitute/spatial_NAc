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
#opt$subset_neurons <- TRUE
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

# Read in the computed t-statistics

if(opt$subset_neurons){
  t_stat_mat_human <-  readRDS(file.path(res_dir, "t_stat_mat_human_neurons.rds"))
}else{
  t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human.rds"))
}

if(opt$subset_neurons){
    t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat_neurons.rds"))
}else{
    t_stat_mat_rat <- readRDS(file.path(res_dir, "t_stat_mat_rat.rds"))
}

# Remove any NA values from human and rat t stats
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

# Plot the heatmap of correlation
colrange <-  seq(-.6,.6, by = 0.01)
colorpal <- colorRampPalette(brewer.pal(n = 7, name = "PRGn"))(length(colrange))

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("_", " ", colnames(corr_mat))

if(opt$subset_neurons){
  row_order <- c("Drd1 MSN 1", "Drd1 MSN 2", "Drd2 MSN 1", "Drd2 MSN 2", "Grm8 MSN", "Glutamatergic", "Pvalb Interneuron", "Chat Interneuron","Sst Interneuron", "GABAergic", "Drd3 MSN")
  col_order <- c("DRD1 MSN A", "DRD1 MSN C", "DRD2 MSN A", "DRD2 MSN B", "DRD1 MSN B", "DRD1 MSN D", "Excitatory", 
  "Inh A", "Inh B", "Inh C", "Inh D", "Inh E", "Inh F")
}else{
  row_order <- c("Drd1 MSN 1", "Drd1 MSN 2", "Drd2 MSN 1", "Drd2 MSN 2", "Grm8 MSN", "Glutamatergic", "Pvalb Interneuron", "Chat Interneuron","Sst Interneuron", "GABAergic", "Drd3 MSN", 
  "Astrocyte", "Microglia", "Polydendrocyte", "Olig 1", "Mural")
  col_order <- c("DRD1 MSN A", "DRD1 MSN C", "DRD2 MSN A", "DRD2 MSN B", "DRD1 MSN B", "DRD1 MSN D", "Excitatory", 
  "Inh A", "Inh B", "Inh C", "Inh D", "Inh E", "Inh F", "Astrocyte A", "Astrocyte B", "Microglia", "OPC", "Oligo", "Endothelial", "Ependymal")
}

corr_mat <- corr_mat[ ,match(col_order, colnames(corr_mat))]
corr_mat <- corr_mat[match(row_order, rownames(corr_mat)), ]

p <- pheatmap(corr_mat,
         color=colorpal,
         cluster_cols=F, 
         cluster_rows=F,
         breaks=colrange,
         fontsize=11, 
         fontsize_row=11, 
         fontsize_col=11,
         fontsize_number=8,
         legend_breaks=c(seq(-.6,.6, by = 0.3)), 
         display_numbers=T, 
         number_format="%.2f", 
         number_color = "black", silent = TRUE)

if(opt$subset_neurons){
  pdf(file.path(plot_dir, "correlation_heatmap_sn_neurons.pdf"), width = 8, height = 5)
  print(p)
  dev.off()
}else{
  pdf(file.path(plot_dir, "correlation_heatmap_sn.pdf"), width = 10, height = 6)
  print(p)
  dev.off()
}
