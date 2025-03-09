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

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "non_human_primate")
plot_dir <- here("plots", "17_cross_species_comparison", "non_human_primate")

subset_neurons <- FALSE
# Read in human snRNA-seq
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
sce <- sce[ ,!sce$CellType.Final == "Neuron_Ambig"]
# Read in the reference data set
ext_dir <- here("processed-data", "17_cross_species_comparison", "NHP_data")

if(subset_neurons){
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_MSNs_processed_final.rds"))
}else{
  sce_ext <- readRDS(file = file.path(ext_dir, "GSE167920_Results_full_nuclei_processed_final.rds"))
}
if(subset_neurons){
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
if(subset_neurons){
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
if(subset_neurons){
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

if(subset_neurons){
    t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human_neurons.rds"))
    t_stat_mat_nhp <- readRDS(file.path(res_dir, "t_stat_mat_nhp_neurons.rds"))

}else{
    t_stat_mat_human <- readRDS(file.path(res_dir, "t_stat_mat_human.rds"))
    t_stat_mat_nhp <- readRDS(file.path(res_dir, "t_stat_mat_nhp.rds"))
}

# Remove any NA values from human and rat t stats
t_stat_mat_nhp <- t_stat_mat_nhp[rowSums(is.na(t_stat_mat_nhp)) == 0, ]
t_stat_mat_human <- t_stat_mat_human[rowSums(is.na(t_stat_mat_human)) == 0, ]

DEGs_df <- DEGs_df[DEGs_df$Gene_symbol %in% rownames(t_stat_mat_human), ]
DEGs_df <- DEGs_df[DEGs_df$Monkey_Symbol %in% rownames(t_stat_mat_nhp), ]

t_stat_mat_nhp <- t_stat_mat_nhp[rownames(t_stat_mat_nhp) %in% DEGs_df$Monkey_Symbol, ]
t_stat_mat_nhp <- t_stat_mat_nhp[match(DEGs_df$Monkey_Symbol, rownames(t_stat_mat_nhp)), ]

t_stat_mat_human <- t_stat_mat_human[rownames(t_stat_mat_human) %in% DEGs_df$Gene_symbol, ]
t_stat_mat_human <- t_stat_mat_human[match(DEGs_df$Gene_symbol, rownames(t_stat_mat_human)), ]

corr_mat <- matrix(NA, nrow = dim(t_stat_mat_nhp)[2],ncol = dim(t_stat_mat_human)[2])
rownames(corr_mat) <- colnames(t_stat_mat_nhp)
colnames(corr_mat) <- colnames(t_stat_mat_human)
for(i in rownames(corr_mat)){
  for(j in colnames(corr_mat)){
    corr_mat[i, j] <- cor(as.numeric(t_stat_mat_nhp[ ,i]), as.numeric(t_stat_mat_human[ ,j]))
  }
}

# Plot the heatmap of correlation
colrange <-  seq(-.6,.6, by = 0.01)
colorpal <- colorRampPalette(brewer.pal(n = 7, name = "PRGn"))(length(colrange))

rownames(corr_mat) <- gsub("\\.", " ", rownames(corr_mat))
colnames(corr_mat) <- gsub("_", " ", colnames(corr_mat))

if(subset_neurons){
  row_order <- c("D1 Matrix", "D2 Matrix", "D1 Striosome", "D2 Striosome", "D1 Shell/OT", "D2 Shell/OT", "D1 ICj", "D1 NUDAP", "D1/D2 Hybrid")
  col_order <- c("DRD1 MSN A", "DRD2 MSN A", "DRD2 MSN B", "DRD1 MSN C", "DRD1 MSN D", "DRD1 MSN B")
}else{
  row_order <- c("DRD1", "DRD2", "Interneurons", "Astrocytes", "Microglia", "Oligos Pre", "Oligos", "Mural/Fibroblast","Endothelial")
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

if(subset_neurons){
  pdf(file.path(plot_dir, "correlation_heatmap_sn_neurons.pdf"), width = 5, height = 4)
  print(p)
  dev.off()
}else{
  pdf(file.path(plot_dir, "correlation_heatmap_sn.pdf"), width = 10, height = 4)
  print(p)
  dev.off()
}
