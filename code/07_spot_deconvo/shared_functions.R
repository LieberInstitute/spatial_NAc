################################################################################
#   Functions related to spot plots
################################################################################

#   Wrapper around spatialLIBD::vis_clus and spatialLIBD::vis_gene, suitable
#   for merged samples (each sample in the SpatialExperiment 'spe' is a donor
#   consisting of multiple capture areas whose spatial coordinates may
#   haphazardly overlap (arbitrary rotations and translations). Designed to
#   address the overplotting and visually awkward spot positions in this case.
#
#   Spot sizes are *almost* consistent among donors, regardless of full-
#   resolution image dimensions, when title is NULL, include_legend is FALSE,
#   and the plot is saved to a square output (e.g. PDF with 7in width and
#   height). However, ggplot does not seem to scale plots of different aspect
#   ratios exactly consistently when writing to PDF (untested for other formats)
#
#   Return a spot plot of sample 'sampleid', assumed to be a donor. 'coldatavar'
#   (character(1)) must be a valid colname in colData(spe).
#
#   Spatial coordinates are fit to a new Visium-like hexagonal grid with the
#   same inter-spot distance (100um in valid Visium data); this prevents
#   overplotting in regions of overlapping capture areas, and ensures all points
#   are equidistant. Where this results in duplicate spot mappings (a pair of
#   spots from two different capture areas round to the same array coordinates),
#   coloring depends on whether the plotted variable is discrete:
#       If [is_discrete], one spot is arbitrarily chosen;
#       If not [is_discrete], the plotted variable is averaged over
#           duplicate-mapped spots and colored accordingly
#
#   spe:            passed to 'spe' in either 'vis_gene' or 'vis_clus'
#   sample_id:      passed to 'sampleid'
#   title:          title for the plot, expected to be one line (avoid use of
#                   "\n")
#   var_name:       passed to 'geneid' for 'vis_gene' and 'clustervar' for
#                   'vis_clus'
#   include_legend: (logical) if FALSE, remove the legend
#   is_discrete:    (logical) if TRUE, use 'vis_clus'; otherwise, use 'vis_gene'
#   colors:         passed to 'colors' for 'vis_gene' if not [is_discrete]
#   assayname:      passed to 'assayname' for 'vis_gene' if not [is_discrete]
#   minCount:       passed to 'minCount' for 'vis_gene' if not [is_discrete]
#
#   Returns a ggplot object
spot_plot <- function(
        spe, sample_id, title, var_name, include_legend, is_discrete,
        colors = NULL, assayname = "logcounts", minCount = 0.5
    ) {
    IDEAL_POINT_SIZE <- 200

    #   Subset to specific sample ID
    spe_small = spe[, spe$sample_id == sample_id]
    
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
    
    #   Find the appropriate spot size for this donor. This can vary because
    #   ggplot downscales a plot the fit desired output dimensions (in this
    #   case presumably a square region on a PDF), and stitched images can vary
    #   in aspect ratio. Also, lowres images always have a larger image
    #   dimension of 1200, no matter how many spots fit in either dimension.
    spot_size = IDEAL_POINT_SIZE * INTER_SPOT_DIST_PX *
       scaleFactors(spe_small) / max(dim(imgData(spe_small)$data[[1]]))

    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   For 'vis_clus' only, supply a color scale if 'color' is not NULL
        if (is.null(colors)) {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = FALSE, point_size = spot_size
            )
        } else {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = FALSE, colors = colors,
                point_size = spot_size
            )
        }
    } else {
        #   For continuous variables, we'll average the values for overlapping
        #   spots
        a = colData(spe_small) |>
            as_tibble() |>
            group_by(array_row, array_col) |>
            mutate(
                !!var_name := mean(!!sym(var_name)),
                is_first_spot = cur_group_id(),
            ) |>
            #   This "extra" mutate is needed due to https://github.com/tidyverse/dplyr/issues/6889
            mutate(is_first_spot = !(duplicated(is_first_spot)))
        
        #   Keep only one spot in cases of multiple spots mapping to the same
        #   array coordinates. Update the colData to use the average over
        #   duplicated spots
        spe_small = spe_small[, a$is_first_spot]
        colData(spe_small)[[var_name]] = a |>
            filter(is_first_spot) |>
            pull({{ var_name }})

        p <- vis_gene(
            spe_small,
            sampleid = sample_id, geneid = var_name, return_plots = TRUE,
            spatial = FALSE, point_size = spot_size, assayname = assayname,
            cont_colors = viridisLite::plasma(21), alpha = 1, auto_crop = FALSE,
            minCount = minCount
        )
    }

    #   Remove the legend if requested
    if (!include_legend) {
        p <- p + theme(legend.position = "none")
    }

    #   Overwrite the title
    p <- p + labs(title = title)

    return(p)
}

#   Set the largest value for the color/fill scale to [upper_limit] and return
#   a copy of the plot 'p' with that modification. 'min_count' should be the
#   value of 'minCount' passed to 'spot_plot', used to create 'p'.
overwrite_scale <- function(p, upper_limit, min_count) {
    p +
        scale_fill_gradientn(
            colors = viridisLite::plasma(21),
            limits = c(0, upper_limit), na.value = c("#CCCCCC40")
        ) +
        scale_color_gradientn(
            colors = viridisLite::plasma(21),
            limits = c(0, upper_limit), na.value = c("#CCCCCC40")
        ) +
        labs(fill = paste(" min >", min_count))
}

#   Given a list of ggplot objects 'plot_list', write up to 3 PDFs with
#   variations of this list. Write the default version as a grid with [ncol]
#   columns, expected to have a legend and a 1-line title; this is designed for
#   best internal viewing. Write a second version with legends but no title, for
#   use in grid-based plots in the manuscript. If [include_individual], write a
#   third version with one plot per page, the existing titles, and no legends.
#
#   plot_list:          list of ggplot objects to plot to multiple PDFs
#   n_col:              (integer) number of columns in grid versions of the
#                       plots
#   plot_dir:           (character) path to base directory for writing plots
#   file_prefix:        (character) start of filename (without extension)
#   include_individual: (logical) if TRUE, write a 3rd version of plots with one
#                       plot per page
#
#   Returns NULL.
write_spot_plots <- function(plot_list, n_col, plot_dir, file_prefix, include_individual) {
    #   Scaling factor for 'plot_grid'. 1.03 seems to reduce whitespace without
    #   introducing overlap, but larger values quickly introduce overlap
    SCALE <- 1.03

    #   Create two alternate versions of the default plot for use in the
    #   manuscript. One version has one plot per page, with consistent
    #   point size and aspect ratio with other individual spot plots. The
    #   other is a grid version with a legend but no title, consistent with
    #   grid plots in the manuscript
    plot_list_individual <- list()
    plot_list_no_title <- list()
    for (i in 1:length(plot_list)) {
        plot_list_individual[[i]] <- plot_list[[i]] +
            theme(legend.position = "none")
        plot_list_no_title[[i]] <- plot_list[[i]] + labs(title = NULL)
    }

    #   Internal-viewing version
    pdf(
        file.path(plot_dir, paste0(file_prefix, ".pdf")),
        width = 7 * n_col,
        height = 7 * length(plot_list) / n_col
    )
    print(plot_grid(plotlist = plot_list, ncol = n_col), scale = SCALE)
    dev.off()

    #   Individual version, if requested
    if (include_individual) {
        pdf(file.path(plot_dir, paste0(file_prefix, "_individual.pdf")))
        print(plot_list_individual)
        dev.off()
    }

    #   No-title version
    pdf(
        file.path(
            plot_dir, paste0(file_prefix, "_no_title.pdf")
        ),
        width = 7 * n_col,
        height = 7 * length(plot_list) / n_col
    )
    print(plot_grid(plotlist = plot_list_no_title, ncol = n_col, scale = SCALE))
    dev.off()
}
