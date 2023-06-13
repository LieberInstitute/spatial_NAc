library(here)
library(tidyverse)
library(jaffelab)

tissue_paths = list.files(
    here('processed-data', '02_image_stitching'),
    pattern = '^tissue_positions_.*(_[ABCD]1){2,4}\\.csv$',
    full.names = TRUE
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

plot_dir = here('plots', '02_image_stitching')

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6
NUM_COL = 128
NUM_ROW = 78

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

################################################################################
#   Overwrite array_row and array_col with new values appropriate for a merged
#   image
################################################################################

#   Essentially, we're placing the merged set of transformed images on a 
#   "new Visium capture area". It contains the traditional number of rows and
#   columns but on a much larger area. By convention, we'll have array row
#   increase with pixel row and array col increase with pixel col.

MIN_ROW = min(these_coords$pxl_row_in_fullres)
INTERVAL_ROW = (max(these_coords$pxl_row_in_fullres) - MIN_ROW) / (NUM_ROW - 1)
    
MIN_COL = min(these_coords$pxl_col_in_fullres)
INTERVAL_COL = (max(these_coords$pxl_col_in_fullres) - MIN_COL) / (NUM_COL - 1)

these_coords$array_row = round(
    (these_coords$pxl_row_in_fullres - MIN_ROW) / INTERVAL_ROW
)
these_coords$array_col = (these_coords$pxl_col_in_fullres - MIN_COL) /
    INTERVAL_COL

#   Actually, columns are constrained to be even for even rows and odd for odd
#   rows. Round accordingly
these_coords$array_col[these_coords$array_row %% 2 == 0] = these_coords |>
    filter(these_coords$array_row %% 2 == 0) |>
    mutate(temp = round(array_col / 2) * 2) |>
    pull(temp)
these_coords$array_col[these_coords$array_row %% 2 == 1] = these_coords |>
    filter(these_coords$array_row %% 2 == 1) |>
    mutate(temp = round(array_col / 2 + 0.5) * 2 - 1) |>
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
stopifnot(max(these_coords$array_row) == 77)
stopifnot(max(these_coords$array_col) == 127)

#   Check an eccentric detail of Visium arrays: (0, 0) cannot exist
if (any((these_coords$array_row == 0) & (these_coords$array_col == 0))) stop()

################################################################################
#   Visually assess array_row and array_col
################################################################################

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

pdf(file.path(plot_dir, paste0('spots_', this_slide, '.pdf')))
print(p)
dev.off()
