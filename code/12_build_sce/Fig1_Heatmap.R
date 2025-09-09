#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(here)

#Load sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

#Remove the neuronal ambiguous population from further analysis. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

#Load cluster colors
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE) 

###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce$CellType.Final)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SYT1","SNAP25",
                 "GAD1","PPP1R1B","BCL11B", #MSN
                 "DRD1","RELN","TAC1","PDYN",#D1-A
                 "RXFP1","CPNE4","OPRM1",#D1-B
                 "TAFA1","RXFP2", #D1-C
                 "VWC2L","NPFFR2", #D1-D
                 "DRD2","ADORA2A","PENK","PTPRM", #D2-A
                 "CLMP","TTN", #D2-B
                 "IL1RAPL2","PDGFD", #Inh_A
                 "VIP","CCK", #Inh_B
                 "GLP1R","TAC3", #Inh_C
                 "SLC5A7","CHAT", #Inh_D
                 "NPY","CORT","SST", #Inh_E
                 "KCNC2","ANK1", #Inh_F
                 "SLC17A7","TBR1", #Exctitatory
                 "AQP4","GFAP","SLC1A2","FGFR3", #Astrocyte_A
                 "TNC","CD44", #Astrocyte_B
                 "FOXJ1","CFAP44", #Ependymal
                 "C3","ARHGAP15",#Microglia
                 "MBP","MOBP", #Oligodendrocyte
                 "PDGFRA","VCAN",#OPC
                 "CLDN5","DCN") #Endothelial
#marker labels
marker_labels <- c(rep("Pan_Neuronal",2),
                   rep("MSN",3),
                   rep("DRD1_MSN_A",4),
                   rep("DRD1_MSN_B",3),
                   rep("DRD1_MSN_C",2),
                   rep("DRD1_MSN_D",2),
                   rep("DRD2_MSN_A",4),
                   rep("DRD2_MSN_B",2),
                   rep("Inh_A",2),
                   rep("Inh_B",2),
                   rep("Inh_C",2),
                   rep("Inh_D",2),
                   rep("Inh_E",3),
                   rep("Inh_F",2),
                   rep("Excitatory",2),
                   rep("Astrocyte_A",4),
                   rep("Astrocyte_B",2),
                   rep("Ependymal",2),
                   rep("Microglia",2),
                   rep("Oligo",2),
                   rep("OPC",2),
                   rep("Endothelial",2))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

cluster_cols <- cluster_cols[-14]

colors_markers <- list(marker = cluster_cols[match(unique(marker_labels),names(cluster_cols))])

colors_markers$marker[1] <- "gray60"
colors_markers$marker[2] <- "gray85"
names(colors_markers$marker)[1] <- "Pan_Neuronal"
names(colors_markers$marker)[2] <- "MSN"

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(DRD1_MSN_A = "DRD1_MSN_A",
                     DRD1_MSN_B = "DRD1_MSN_B",
                     DRD1_MSN_C = "DRD1_MSN_C",
                     DRD1_MSN_D = "DRD1_MSN_D",
                     DRD2_MSN_A = "DRD2_MSN_A",
                     DRD2_MSN_B = "DRD2_MSN_B",
                     Inh_A = "Inh_A",
                     Inh_B = "Inh_B",
                     Inh_C = "Inh_C",
                     Inh_D = "Inh_D",
                     Inh_E = "Inh_E",
                     Inh_F = "Inh_F",
                     Excitatory = "Excitatory",
                     Astrocyte_A = "Astrocyte_A",
                     Astrocyte_B = "Astrocyte_B",
                     Ependymal = "Ependymal",
                     Microglia = "Microglia",
                     Oligo = "Oligo",
                     OPC = "OPC",
                     Endothelial = "Endothelial")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = cluster_cols[names(cluster_pops)])

row_ha <- rowAnnotation(population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = FALSE,
                        col = pop_markers)



dat <- assay(sce,"logcounts")
rownames(dat) <- rowData(sce)$gene_name
dim(dat)
#[1]  36601 103339

dat <- dat[markers_all,]
dim(dat)
#[1]     51 103339

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]

max(hm_mat)
#[1] 4.248515

min(hm_mat)
#[1] -1.478196

col_fun <- circlize::colorRamp2(c(-1.5,0,4.5),c("blue","white","red"))

Highlight_Genes <- c("SYT1","SNAP25",
                     "GAD1","PPP1R1B","BCL11B", #MSN
                     "DRD1","PDYN",
                     "RXFP1","CPNE4","OPRM1", #D1-B
                     "TAFA1","RXFP2", #D1-C
                     "VWC2L","NPFFR2", #D1-D
                     "DRD2","ADORA2A","PENK",
                     "VIP","CCK", #Inh_B
                     "GLP1R","TAC3", #Inh_C
                     "SLC5A7","CHAT", #Inh_D
                     "NPY","CORT","SST", #Inh_E
                     "SLC17A7","TBR1", #Exctitatory
                     "AQP4","GFAP","SLC1A2", #Astrocyte_A
                     "FOXJ1","CFAP44", #Ependymal
                     "C3","ARHGAP15",#Microglia
                     "MBP","MOBP", #Oligodendrocyte
                     "PDGFRA","VCAN",#OPC
                     "CLDN5","DCN") #Endothelial

# Set all row names to black by default
gene_colors <- rep("black", length(markers_all))
names(gene_colors) <- markers_all

# Set specific genes to red
gene_colors[Highlight_Genes] <- "red"

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
                              column_names_gp = gpar(col = gene_colors),
                              rect_gp = gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)


pdf(file = here("plots","12_snRNA","Expression","CellType_Final_Heatmap_BroadMarkers.pdf"),
    width = 16,
    height = 12)
draw(hm)
dev.off()


sessionInfo()
