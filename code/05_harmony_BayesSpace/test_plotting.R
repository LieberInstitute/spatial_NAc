library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe.rds')
sample_id = "Br2720"
cluster_var = "X10x_kmeans_7_clusters"
plot_path = here('plots', '05_harmony_BayesSpace', 'test.pdf')

spe = readRDS(spe_path)

# colnames(spe) = spe$key
# spe$pxl_row_in_fullres_original <- spe$pxl_row_in_fullres
# spe$pxl_col_in_fullres_original <- spe$pxl_col_in_fullres
# spe$pxl_row_in_fullres <- spe$pxl_col_in_fullres <- NULL

vis_clus_merged = function(spe, sampleid, clustervar) {
    #   Subset to specific sample ID
    spe_small = spe[, spe$sample_id == sampleid]
    
    #   Determine some pixel values for the horizontal and vertical bounds
    #   of the spots
    MIN_ROW = min(spatialCoords(spe_small)[, 'pxl_col_in_fullres'])
    MIN_COL = min(spatialCoords(spe_small)[, 'pxl_row_in_fullres'])
    MAX_COL = max(spatialCoords(spe_small)[, 'pxl_row_in_fullres'])

    #   The distance between spots (in pixels) is double the average distance
    #   between array columns 
    INTER_SPOT_DIST_PX = 2 * (MAX_COL - MIN_COL) /
        (max(spe_small$array_col) - min(spe_small$array_col))
    
    #   Find the distance (in pixels) between array columns and between array
    #   rows, respectively
    INTERVAL_COL = INTER_SPOT_DIST_PX / 2
    INTERVAL_ROW = INTER_SPOT_DIST_PX * cos(pi / 6)
    
    #   Now overwrite spatialCoords to use pixel coordinates rounded to the
    #   nearest array value. This ensures smooth and equidistant plotted points
    #   even in overlapping capture areas
    spatialCoords(spe_small)[, 'pxl_col_in_fullres'] = MIN_ROW +
        spe_small$array_row * INTERVAL_ROW
    spatialCoords(spe_small)[, 'pxl_row_in_fullres'] = MAX_COL -
        spe_small$array_col * INTERVAL_COL

    #   Return a plot from 'vis_clus' that uses the rounded spatialCoords
    p = vis_clus(
        spe_small, sampleid = sampleid, clustervar = clustervar,
        auto_crop = FALSE, point_size = 1
    )
    return(p)
}

p = vis_clus_merged(spe, sampleid = sample_id, clustervar = cluster_var)
pdf(plot_path)
print(p)
dev.off()
