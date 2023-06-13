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

#   55-micrometer diameter for Visium spot
SPOT_DIAMETER_M = 55e-6
NUM_COL = 128
NUM_ROW = 78

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
this_sample_info = sample_info |>
    filter(slide_num == 'V12D07-074')

these_coords = coords |> filter(sample_id %in% this_sample_info$sample_id)

################################################################################
#   Overwrite array_row and array_col with new values appropriate for a merged
#   image
################################################################################

#   Essentially, we're placing the merged set of transformed images on a 
#   "new Visium capture area". It contains the traditional number of rows and
#   columns but on a much larger area. By convention, we'll have array row
#   increase with pixel row and array col increase with pixel col.

#-------------------------------------------------------------------------------
#   Determine which array row a given spot should be assigned
#-------------------------------------------------------------------------------

MAX = max(these_coords$pxl_row_in_fullres)
MIN = min(these_coords$pxl_row_in_fullres)
these_coords$array_row = round(
    (NUM_ROW - 1) * (these_coords$pxl_row_in_fullres - MIN) / (MAX - MIN)
)

#-------------------------------------------------------------------------------
#   Given the newly assigned array row, determine which array column a given
#   spot should be assigned
#-------------------------------------------------------------------------------

MAX = max(these_coords$pxl_col_in_fullres)
MIN = min(these_coords$pxl_col_in_fullres)

#   The pixel distance between array row i and row i+1
INTERVAL_SIZE = (MAX - MIN) / (NUM_COL - 1)

#   'min_vec' gives the minimum 'pxl_row_in_fullres' values possible for spots
#   in the corresponding new array row. Odd rows have columns shifted in the
#   positive direction by one INTERVAL_SIZE. Also, row 0 (strangely) essentially
#   skips the first column (this is just how Visium capture areas are)
min_vec = MIN + INTERVAL_SIZE * (these_coords$array_row %% 2 == 1)
min_vec[match(0, these_coords$array_row)] = MIN + INTERVAL_SIZE * 2

#   The interval size here is doubled because for a given row, only every other
#   column value is possible. We multiply by 2 because we're still indexing as
#   if every column is possible
these_coords$array_col = 2 * round(
        (these_coords$pxl_col_in_fullres - min_vec) / (2 * INTERVAL_SIZE)
    ) +
    #   Odd rows have columns start at 1, not 0
    (these_coords$array_row %% 2 == 1)

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
