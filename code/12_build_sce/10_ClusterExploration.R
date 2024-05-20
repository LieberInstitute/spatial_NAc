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
