library(here)
library(tidyverse)
library(jaffelab)
library(rjson)
library(cowplot)

tissue_paths = list.files(
    here('processed-data', '02_image_stitching'),
    pattern = '^tissue_positions_.*(_[ABCD]1){2,4}\\.csv$',
    full.names = TRUE
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

plot_dir = here('plots', '02_image_stitching')

#   55-micrometer diameter for Visium spot; 100 micrometers between spots;
#   65-micrometer spot diameter used in 'spot_diameter_fullres' calculation for
#   spaceranger JSON. See documentation for respective quantities. The
#   difference between 55 and 65 does indeed exist and is properly documented,
#   but is likely a bug in the sense the choice was probably unintentional
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487812-What-is-the-size-of-the-spots-on-the-Visium-Gene-Expression-Slide-
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487892-How-much-space-is-there-between-spots-referred-to-as-white-space-
#   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
SPOT_DIAMETER_M = 55e-6
SPOT_DIAMETER_JSON_M = 65e-6
INTER_SPOT_DIST_M = 100e-6

################################################################################
#   Functions
################################################################################

#   Given 'these_coords', a tibble whose rows represent samples of the same
#   donor, and containing columns 'array_row', 'array_col',
#   'pxl_row_in_fullres', and 'pxl_col_in_fullres', modify the 'array_row' and
#   'array_col' columns to represent a larger Visium capture area containing
#   all samples in a common coordinate system. 'inter_spot_dist_px' gives the
#   distance between any 2 spots in the new coordinates, and so the number of
#   array rows/cols generally changes from the Visium standards of 78 and 128
#   (and even changes in ratio between num rows and num cols). Return the same
#   tibble with modified 'array_row' + 'array_col' columns, as well as new
#   'pxl_row_rounded' and 'pxl_col_rounded' columns representing the pixel
#   coordinates rounded to the nearest exact array coordinates.
fit_to_array = function(these_coords, inter_spot_dist_px) {
    MIN_ROW = min(these_coords$pxl_row_in_fullres)
    INTERVAL_ROW = inter_spot_dist_px * cos(pi / 6)
        
    MIN_COL = min(these_coords$pxl_col_in_fullres)
    INTERVAL_COL = inter_spot_dist_px / 2

    these_coords$array_row = round(
        (these_coords$pxl_row_in_fullres - MIN_ROW) / INTERVAL_ROW
    )
    these_coords$array_col = (these_coords$pxl_col_in_fullres - MIN_COL) /
        INTERVAL_COL

    #   Actually, columns are constrained to be even for even rows and odd for
    #   odd rows. Round accordingly
    these_coords$array_col[these_coords$array_row %% 2 == 0] = these_coords |>
        filter(these_coords$array_row %% 2 == 0) |>
        mutate(temp = round(array_col / 2) * 2) |>
        pull(temp)
    these_coords$array_col[these_coords$array_row %% 2 == 1] = these_coords |>
        filter(these_coords$array_row %% 2 == 1) |>
        mutate(temp = round(array_col / 2 - 0.5) * 2 + 1) |>
        pull(temp)

    #   Now make new pixel columns based on just the array values (these columns
    #   give the coordinates for given array row/cols)
    these_coords$pxl_row_rounded = MIN_ROW + these_coords$array_row * INTERVAL_ROW
    these_coords$pxl_col_rounded = MIN_COL + these_coords$array_col * INTERVAL_COL

    #-------------------------------------------------------------------------------
    #   Verify the newly assigned array row and cols have reasonable values
    #-------------------------------------------------------------------------------

    #   Even array rows can only use even column indices
    these_coords |>
        filter(array_row %% 2 == 0) |>
        summarize(a = all(array_col %% 2 == 0)) |>
        pull(a) |>
        stopifnot()

    #   Odd array rows can only use odd column indices
    these_coords |>
        filter(array_row %% 2 == 1) |>
        summarize(a = all(array_col %% 2 == 1)) |>
        pull(a) |>
        stopifnot()

    #   Check range of array row and col
    stopifnot(min(these_coords$array_row) == 0)
    stopifnot(min(these_coords$array_col) == 0)

    #   Check an eccentric detail of Visium arrays: (0, 0) cannot exist
    if (any((these_coords$array_row == 0) & (these_coords$array_col == 0))) stop()

    return(these_coords)
}

################################################################################
#   Read in sample info and spot coordinate info
################################################################################

coords = do.call(rbind, lapply(tissue_paths, read.csv)) |> as_tibble()
coords$sample_id = paste(
    ss(coords$key, '_', 2), ss(coords$key, '_', 3), sep = '_'
)

sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(
        sample_id = X,
        array_num = 'Array..',
        slide_num = 'Slide..' 
    )

print("Correlation relationship between array (a) and pixel (p) row/col:")
coords |>
    group_by(sample_id) |>
    summarize(
        arow_v_prow = cor(array_row, pxl_row_in_fullres),
        arow_v_pcol = cor(array_row, pxl_col_in_fullres),
        acol_v_prow = cor(array_col, pxl_row_in_fullres),
        acol_v_pcol = cor(array_col, pxl_col_in_fullres)
    ) |>
    print()

#   Subset to 1 donor (TODO: not by slide; make more general)
this_slide = 'V12D07-074'
this_sample_info = sample_info |>
    filter(slide_num == this_slide)

these_coords = coords |> filter(sample_id %in% this_sample_info$sample_id)

#-------------------------------------------------------------------------------
#   Select an appropriate number of array rows and columns for this slide such
#   that the distance between spots on the new grid matches the Visium standard
#-------------------------------------------------------------------------------

sr_json = fromJSON(
    file = file.path(
        this_sample_info$spaceranger_dir[1], "scalefactors_json.json"
    )
)

#   Since we have the spot diameter both in pixels and meters, compute the
#   image's pixel/m ratio. Then use that to compute the distance between spots
#   in pixels
PX_PER_M = sr_json$spot_diameter_fullres / SPOT_DIAMETER_JSON_M
INTER_SPOT_DIST_PX = INTER_SPOT_DIST_M * PX_PER_M

#   Sanity check on existing coordinates: spots should be 100um apart. If not,
#   modify. TODO: look for pairs of array coordinates that exist; don't hardcode
#   10, 11
a = these_coords |>
    filter(array_row == 10, array_col == 80) |>
    head(1) |>
    select(pxl_row_in_fullres, pxl_col_in_fullres)

b = these_coords |>
    filter(array_row == 11, array_col == 81) |>
    head(1) |>
    select(pxl_row_in_fullres, pxl_col_in_fullres)

tol = 2
observed_dist = sqrt(rowSums((a - b) ** 2))
if (abs(observed_dist - INTER_SPOT_DIST_PX) > tol) {
    stop('Transformed spots output from ImageJ/Samui were not close to 100um apart!')
}

################################################################################
#   Visually assess array_row and array_col after fitting to "new grid"
################################################################################

#   Adjust 'array_row' and 'array_col' with values appropriate for the new
#   coordinate system (a larger Visium grid with equal inter-spot distances)
these_coords = fit_to_array(these_coords, INTER_SPOT_DIST_PX)

#-------------------------------------------------------------------------------
#   Spot plots
#-------------------------------------------------------------------------------

#   Plot the transformed spots as-is for reference
p = ggplot(these_coords) +
    geom_point(
        aes(
            x = pxl_col_in_fullres,
            y = max(pxl_row_in_fullres) - pxl_row_in_fullres,
            color = sample_id
        ),
        size = 0.1
    ) +
    coord_fixed() +
    labs(
        x = 'Full-res pixels (width)', y = 'Full-res pixels (height)',
        color = 'Sample ID'
    ) + 
    guides(color = guide_legend(override.aes = list(size = 1.5)))

pdf(file.path(plot_dir, paste0('raw_spots_', this_slide, '.pdf')))
print(p)
dev.off()

#   Plot the spots aligned to the new Visium grid
p = ggplot(these_coords) +
    geom_point(
        aes(
            x = pxl_col_rounded,
            y = max(pxl_row_rounded) - pxl_row_rounded,
            color = sample_id
        ),
        size = 0.1
    ) +
    coord_fixed() +
    labs(
        x = 'Full-res pixels (width)', y = 'Full-res pixels (height)',
        color = 'Sample ID', title = 'Spot Positions (Aligned to Array)'
    ) + 
    guides(color = guide_legend(override.aes = list(size = 1.5)))

pdf(file.path(plot_dir, paste0('aligned_spots_', this_slide, '.pdf')))
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   Plots measuring mapping of multiple spots to the same array coordinates
#-------------------------------------------------------------------------------

#   Verify that for a single sample, mappings were unique
dup_coords = these_coords |>
    group_by(sample_id, array_col, array_row) |>
    summarize(n = n()) |>
    ungroup()

p1 = ggplot(dup_coords) +
    geom_histogram(
        aes(x = n, fill = sample_id), bins = length(unique(dup_coords$n)),
        color = "black"
    ) +
    facet_wrap(~sample_id, nrow = 1) +
    labs(
        x = "Num spots per array coordinate", fill = "Sample ID",
        title = "Within-sample duplication of spot mappings"
    )

#   Then check how many spots mapped to the same array coordinates, agnostic
#   to sample
dup_coords = these_coords |>
    group_by(array_col, array_row) |>
    summarize(n = n()) |>
    ungroup()

p2 = ggplot(dup_coords) +
    geom_histogram(aes(x = n), bins = length(unique(dup_coords$n))) +
    labs(
        x = "Num spots per array coordinate",
        title = "Overall duplication of spot mappings"
    )

pdf(file.path(plot_dir, 'spot_duplication.pdf'))
plot_grid(p1, p2, nrow = 2)
dev.off()

################################################################################
#   Measure error in aligning (rounding) pixel coordinates to fit array
################################################################################

#   Mean alignment error in spot diameters
print("Mean alignment error in spot diameters:")
these_coords |>
    summarize(
        d = mean(
            (pxl_row_in_fullres - pxl_row_rounded) ** 2 +
            (pxl_col_in_fullres - pxl_col_rounded) ** 2
        ) ** 0.5 / sr_json$spot_diameter_fullres
    ) |>
    pull(d) |>
    round(2) |>
    print()
