#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(33): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(5): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# umap_HARMONY
# mainExpName: NULL
# altExpNames(0):

sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                        "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                                        "Inh_C","Inh_D","Inh_E","Inh_F",
                                        "Inh_G","Inh_H","Inh_I","Excitatory",
                                        "Astrocyte_A","Astrocyte_B","Ependymal","Oligo",
                                        "OPC","Microglia","Macrophage","T-Cell",
                                        "Mural_Endothelial_Fibroblast"))

sce$CellType.Final <- factor(x = sce$CellType.Final,
                             levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                        "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                                        "Inh_C","Inh_D","Inh_E","Inh_F",
                                        "Inh_G","Inh_H","Inh_I","Excitatory",
                                        "Astrocyte_A","Astrocyte_B","Ependymal","Oligo",
                                        "OPC","Microglia","Macrophage","T-Cell",
                                        "Mural_Endothelial_Fibroblast"))


####Remake BrainID by cluster with new order. 
#Calculate cluster percentage by BrainID
brain_by_cluster <- as.data.frame.matrix(table(sce$Brain_ID,sce$CellType.Final))

#Calculate percentages 
brain_by_cluster_pct <- sweep(brain_by_cluster,MARGIN = 2,colSums(brain_by_cluster),"/") * 100
brain_by_cluster_pct$brain <- rownames(brain_by_cluster_pct)

#Melt dataframe and plot
brain_by_cluster_pct_melt <- reshape2::melt(brain_by_cluster_pct)

#make some new brain colors
brain_cols <- Polychrome::createPalette(length(unique(sce$Brain_ID)),
                                        c("#D81B60", "#1E88E5","#004D40"))
names(brain_cols) <- unique(sce$Brain_ID)

brain_cluster_bar <- ggplot(data = brain_by_cluster_pct_melt,aes(x = variable,y = value,fill = brain)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = brain_cols) +
  labs(x = "CellType",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = brain_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","BrainID_by_CellType_Final_Ordered.pdf"),
       height = 12,
       width = 12)

#Remake the above plot with donor on the x-axis. 
cluster_by_brain <- as.data.frame.matrix(table(sce$CellType.Final,sce$Brain_ID))

#Calculate percentages 
cluster_by_brain_pct <- sweep(cluster_by_brain,MARGIN = 2,colSums(cluster_by_brain),"/") * 100
cluster_by_brain_pct$CellType <- rownames(cluster_by_brain_pct)

#Melt dataframe and plot
cluster_by_brain_pct_melt <- reshape2::melt(cluster_by_brain_pct)

#load the colors
load(here("processed-data","12_snRNA","CellType_Final_Cluster_Cols.rda"),verbose = TRUE)
# Loading objects:
#   cluster_cols

cluster_brain_bar <- ggplot(data = cluster_by_brain_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Donor",
       y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(plot = cluster_brain_bar,
       filename = here("plots","12_snRNA","Dim_Red","CellType_Final_By_BrainID_Ordered.pdf"),
       height = 12,
       width = 12)


####Sort by celltype
#Calculate cluster percentage by sort
sort_by_cluster <- as.data.frame.matrix(table(sce$Sort,sce$CellType.Final))

#Calculate percentages 
sort_by_cluster_pct <- sweep(sort_by_cluster,MARGIN = 2,colSums(sort_by_cluster),"/") * 100
sort_by_cluster_pct$Sort <- rownames(sort_by_cluster_pct)

#Melt dataframe and plot
sort_by_cluster_pct_melt <- reshape2::melt(sort_by_cluster_pct)

#make some new sort colors
sort_cols <- Polychrome::createPalette(length(unique(sce$Sort)),
                                        c("#D81B60","#004D40"))
names(sort_cols) <- unique(sce$Sort)

sort_cluster_bar <- ggplot(data = sort_by_cluster_pct_melt,aes(x = variable,y = value,fill = Sort)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = sort_cols) +
  labs(x = "CellType",
       y = "Percent") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = sort_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","Sort_by_CellType_Final_Ordered.pdf"),
       height = 12,
       width = 12)

#Remake above plot but with sort on the x-axis and celltype proportion on y-axis
cluster_by_sort <- as.data.frame.matrix(table(sce$CellType.Final,sce$Sort))

#Calculate percentages 
cluster_by_sort_pct <- sweep(cluster_by_sort,MARGIN = 2,colSums(cluster_by_sort),"/") * 100
cluster_by_sort_pct$CellType <- rownames(cluster_by_sort_pct)

#Melt dataframe and plot
cluster_by_sort_pct_melt <- reshape2::melt(cluster_by_sort_pct)

#Make the bargraph.
sort_cluster_bar <- ggplot(data = cluster_by_sort_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Sort",
       y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(plot = sort_cluster_bar,
       filename = here("plots","12_snRNA","Dim_Red","CellType_By_Sort_Final_Ordered.pdf"),
       height = 12,
       width = 12)


###Save umap as PDF
umap_harmony_pdf <- plotReducedDim(object = sce,dimred = "umap_HARMONY",
               colour_by = "CellType.Final") +
  scale_color_manual(values = cluster_cols) +
  theme_void() +
  theme(legend.position = "none") 

ggsave(plot = umap_harmony_pdf,
       filename = here("plots","12_snRNA","Dim_Red","umap_HARMONY_Final_noText.pdf"),
       height = 12,
       width = 12)

#Also as png. 
ggsave(plot = umap_harmony_pdf,
       filename = here("plots","12_snRNA","Dim_Red","umap_HARMONY_Final_noText.png"),
       height = 12,
       width = 12)

#Violin of general marker genes. 
gen_vln <- plotExpression(object = sce,
               features = c("RBFOX3","GJA1","MOBP","PDGFRA","C3","CD163","GRAP2","EBF1"),
               swap_rownames = "gene_name",
               x = "CellType.Final", colour_by = "CellType.Final",
               ncol = 1,) +
  scale_color_manual(values = cluster_cols) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  theme_void() +
  theme(legend.position = "none") 

ggsave(plot = gen_vln,
       filename = here("plots","12_snRNA","Expression","General_Violin_Figure1.pdf"))

ggsave(plot = gen_vln,
       filename = here("plots","12_snRNA","Expression","General_Violin_Figure1.png"))


### Read in the DEGs
load(here("processed-data","markers_1vAll_CellType_Final.rda"),verbose = TRUE)
# Loading objects:
#   markers_1vALL_df

###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$CellType.Final)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("RBFOX3","SNAP25","SYT1", #Pan neuronal  3
                 "GAD1","GAD2","SLC32A1", #GABAergic
                 "PPP1R1B","BCL11B","RARB", #MSN
                 "DRD1","DRD2", #Dopamine receptor MSNs
                 "SLC17A7", #Exctitatory
                 "AQP4","GFAP", #Astrocyte
                 "FOXJ1","CFAP44", #Ependymal
                 "C3","ARHGAP15",#Microglia
                 "MBP","MOBP", #Oligodendrocyte
                 "PDGFRA",#OPC
                 "CD163", #Macrophage
                 "CD2", #T-Cell
                 "CLDN5") #Mural_Endothelial_Fibroblast

#marker labels
marker_labels <- c(rep("Neuronal",3),
                   rep("GABAergic",3),
                   rep("MSN",5),
                   rep("Excitatory",1),
                   rep("Astrocyte",2),
                   rep("Ependymal",2),
                   rep("Microglia",2),
                   rep("Oligo",2),
                   rep("OPC",1),
                   rep("Other",3))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c(Neuronal = "black",
                                  GABAergic = "#8C564B",
                                  MSN = "#1F77B4",
                                  Excitatory = "#D81B60",
                                  Astrocyte = "#FFC107",
                                  Ependymal = "#009E73",
                                  Microglia = "#88CCEE",
                                  Oligo = "#2af7db",
                                  OPC = "#ff9c1e",
                                  Other = "#b11eff"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = FALSE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(GABAergic = c("Inh_A","Inh_B","Inh_C",
                                   "Inh_D","Inh_E","Inh_F",
                                   "Inh_G","Inh_H","Inh_I"),
                     MSN = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C",
                             "DRD1_MSN_D","DRD2_MSN_A","DRD2_MSN_B"),
                     Excitatory = "Excitatory",
                     Astrocyte = c("Astrocyte_A","Astrocyte_B"),
                     Ependymal = "Ependymal",
                     Microglia = "Microglia",
                     Oligo = "Oligo",
                     OPC = "OPC",
                     Other = c("Macrophage","T-Cell",
                               "Mural_Endothelial_Fibroblast"))

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
#cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(cluster_pops_rev %in% c("GABAergic","MSN","Excitatory"),
                      "Neuron",
                      "Non-neuron")

neuron_pops <- factor(x = neuron_pops,levels = c("Neuron","Non-neuron"))

colors_neurons <- list(class = c(Neuron = "black",
                                 `Non-neuron` = "gray65"))

n <- table(sce$CellType.Final)

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = c(GABAergic = "#8C564B",
                                   MSN = "#1F77B4",
                                   Excitatory = "#D81B60",
                                   Astrocyte = "#FFC107",
                                   Ependymal = "#009E73",
                                   Microglia = "#88CCEE",
                                   Oligo = "#2af7db",
                                   OPC = "#ff9c1e",
                                   Other = "#b11eff"))


row_ha <- rowAnnotation(n = anno_barplot(as.numeric(n[names(cluster_pops_rev)]), 
                                         gp = gpar(fill = "navy"), 
                                         border = FALSE),
                        class = neuron_pops,
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        col = c(pop_markers,colors_neurons))


dat <- assay(sce,"logcounts")
rownames(dat) <- rowData(sce)$gene_name
dim(dat)

dat <- dat[markers_all,]
dim(dat)
dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))))
hm_mat <- hm_mat[names(cluster_pops_rev),]

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5))


png(filename = here("plots","12_snRNA","Expression","CellType_Final_Heatmap_BroadMarkers.png"),
    width = 12,
    height = 12,units = "in",res = 100)
draw(hm)
dev.off()
