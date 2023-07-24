#   Wrapper around spatialLIBD::vis_clus and spatialLIBD::vis_gene, suitable
#   for merged samples (each sample in the SpatialExperiment 'spe' is a donor
#   consisting of multiple capture areas whose spatial coordinates may
#   haphazardly overlap (arbitrary rotations and translations). Designed to
#   address the overplotting and visually awkward spot positions in this case.
#
#   Return a spot plot of sample 'sampleid', assumed to be a donor. 'coldatavar'
#   (character(1)) must be a valid colname in colData(spe). If it is a factor,
#   a discrete color scale is used by calling 'vis_clus'; otherwise, it's
#   assumed to be continuous and calls 'vis_gene'.
#
#   Spatial coordinates are fit to a new Visium-like hexagonal grid with the
#   same inter-spot distance (100um in valid Visium data); this prevents
#   overplotting in regions of overlapping capture areas, and ensures all points
#   are equidistant. Where this results in duplicate spot mappings (a pair of
#   spots from two different capture areas round to the same array coordinates),
#   coloring depends on whether the plotted variable is discrete:
#       If discrete, one spot is arbitrarily chosen;
#       If continuous, the plotted variable is averaged over duplicate-mapped
#           spots and colored accordingly
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
            group_by(array_row, array_col) |>
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
