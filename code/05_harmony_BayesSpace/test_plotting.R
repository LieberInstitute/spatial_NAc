library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe.rds')
plot_dir = here('plots', '05_harmony_BayesSpace')

sample_id = "Br2720"
discrete_var = "X10x_kmeans_7_clusters"
cont_var = "sum_gene"

spe = readRDS(spe_path)

# colnames(spe) = spe$key
# spe$pxl_row_in_fullres_original <- spe$pxl_row_in_fullres
# spe$pxl_col_in_fullres_original <- spe$pxl_col_in_fullres
# spe$pxl_row_in_fullres <- spe$pxl_col_in_fullres <- NULL

vis_merged = function(spe, sampleid, coldatavar) {
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
    
    if (is.factor(spe_small[[coldatavar]])) {
        #   Return a plot from 'vis_clus' that uses the rounded spatialCoords
        p = vis_clus(
            spe_small, sampleid = sampleid, clustervar = coldatavar,
            auto_crop = FALSE, point_size = 1
        )
        return(p)
    } else {
        #   For continuous variables, we'll average the values for overlapping
        #   spots
        a = colData(spe_small) |>
            as_tibble() |>
            group_by(array_num, array_row, array_col) |>
            mutate(
                !!coldatavar := mean(!!sym(coldatavar)),
                is_first_spot = cur_group_id(),
            ) |>
            #   This "extra" mutate is needed due to https://github.com/tidyverse/dplyr/issues/6889
            mutate(is_first_spot = !(duplicated(is_first_spot)))
        
        #   Keep only one spot in cases of multiple spots mapping to the same
        #   array coordinates. Update the colData to use the average over
        #   duplicated spots
        spe_small = spe_small[, a$is_first_spot]
        colData(spe_small)[[coldatavar]] = a |>
            filter(is_first_spot) |>
            pull({{ coldatavar }})
        
        p = vis_gene(
            spe_small, sampleid = sampleid, geneid = coldatavar,
            auto_crop = FALSE, point_size = 1
        )
        return(p)
    } 
}

p = vis_merged(spe, sampleid = sample_id, coldatavar = discrete_var)
pdf(file.path(plot_dir, 'vis_merged_discrete.pdf'))
print(p)
dev.off()

p = vis_merged(spe, sampleid = sample_id, coldatavar = cont_var)
pdf(file.path(plot_dir, 'vis_merged_continuous.pdf'))
print(p)
dev.off()
