# Include required libraries
library(here)
library(SingleCellExperiment)
library(jaffelab)
library(scater)
library(scran)
library(readxl)
library(Polychrome)
library(cluster)
library(limma)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(ggrastr)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(Polychrome)
library(RColorBrewer)
library(forcats)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
data(palette36)
library(getopt)

# Read in the SPE data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_dimRed_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_dir)

# Add final cluster annotations
saveDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
final_clusters <- read.csv(file.path(saveDir, "precast_clusters.csv"))

spe <- spe[ ,colnames(spe) %in% final_clusters$key]
final_clusters <- final_clusters[match(colnames(spe), final_clusters$key), ]
spe$spatialDomains <- final_clusters$cluster
spe$spatialDomains <- factor(spe$spatialDomains, 
levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM"))

###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$spatialDomains)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("GABRQ", "RXFP1", "FOXP2", "SEMA5B","OPRM1", #D1 islands
                 "NGFR","ACTA2", "CLDN5", "SUSD2",  # Endothelial/ Ependymal
                 "SLC17A7", "TBR1", "NEUROD2", "NRN1", #Excitatory
                 "NPY", "CORT", "SST", "KIT", "CHODL", "SLC5A7","GRIN2D","OPRD1","CHRNA3", # Inhibitory
                 "MCTP2", "ADORA2A", "ANO3", "PDE10A", "KCNH4","CALB1","PPP1R1A","PPP1R1B", "TAC1", 
                 "ARHGAP36", "CARTPT", "CALCR", "PNMA5", "NPY2R", "CBLN4", "KCTD4","OPRK1","CNR1",
                 "PDYN", "PENK", "DRD1", "DRD2",  #MSN
                 "MBP", "MOBP", "GJB1", "PLP1", "OPALIN" # White matter
                 )


#marker labels
marker_labels <- c(rep("D1 islands",5),
                   rep("Endothelial/Ependymal",4),
                   rep("Excitatory",4),
                   rep("Inhibitory",9),
                   rep("MSN",22),
                   rep("WM",5))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c("D1 islands" = "#E7298A",
                                  "Endothelial/Ependymal" = "#A6761D",
                                  "Excitatory" = "#D95F02",
                                  "Inhibitory" = "#E6AB02",
                                  "MSN" = "#3288bd",
                                  "WM" = "#666666"))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers, 
                                           annotation_legend_param = list(marker = list(title = "Marker gene class", labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14))))

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(D1_islands = "D1 islands",
                     Endo = c("Endothelial/Ependymal"),
                     Excitatory = "Excitatory",
                     Inhibitory="Inhibitory", 
                     MSN_1 = "MSN 1", 
                     MSN_2 = "MSN 2", 
                     MSN_3 = "MSN 3", 
                     WM = "WM")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
#cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))


#row annotation dataframe. 
# row annotation
pop_markers <- list(population = c(D1_islands = "#E7298A",
                                   Endo = "#A6761D",
                                   Excitatory = "#D95F02",
                                   Inhibitory = "#E6AB02",
                                   MSN_1 = "#66A61E",
                                   MSN_2 = "#1B9E77",
                                   MSN_3 = "#7570B3",
                                   WM = "#666666"))

row_ha <- rowAnnotation(population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = FALSE,
                        col = pop_markers, 
                        annotation_legend_param = list(population = list(ncol = 1,title = "Spatial domain", title_position = "topcenter", labels_gp = gpar(fontsize = 12, fontface = 'bold'), title_gp = gpar(fontsize = 16, fontface = 'bold'))))


dat <- assay(spe,"logcounts")
rownames(dat) <- rowData(spe)$gene_name
dim(dat)
#[1]  36601 103339

dat <- dat[rownames(dat) %in% markers_all,]
dim(dat)
#[1]     42 103339

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]
hm_mat <- hm_mat[ ,match(markers_all, colnames(hm_mat))]

col_fun <- circlize::colorRamp2(c(-2,0,4),c("blue","white","red"))



hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "Centered,scaled expression",
                              column_title = NULL,
                              column_names_gp = gpar(fontface = "italic"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 1),
                              col = col_fun, 
                              border = TRUE,
                              heatmap_legend_param = list(legend_direction = "horizontal",legend_width = unit(6, "cm"), title_position = "topcenter", title_gp = gpar(fontsize = 14), border = "black"))


pdf(file = here("plots","10_post_clustering_analysis","01_pseudobulk_markers","SpatialDomains_Final_Heatmap_BroadMarkers.pdf"),
    width = 15,
    height = 5)
draw(hm, heatmap_legend_side="bottom")
dev.off()