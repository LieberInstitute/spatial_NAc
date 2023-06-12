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

coords = do.call(read.csv, as.list(tissue_paths)) |> as_tibble()
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

#   Subset to 1 donor (TODO: not by slide; make more general)
this_sample_info = sample_info |>
    filter(slide_num == 'V12D07-074')

these_coords = coords |> filter(sample_id %in% this_sample_info$sample_id)

#   TODO: (probably in another script) check if the relationship between
#   array_row and pixel_row, etc is stable. Try multiple datasets and samples

#   Horribly confusingly, 'array_row' increases with 'pxl_col_in_fullres', not
#   'pxl_row_in_fullres'
MAX = max(these_coords$pxl_col_in_fullres)
MIN = min(these_coords$pxl_col_in_fullres)
these_coords$array_row_new = round(
    (NUM_ROW - 1) * (these_coords$pxl_col_in_fullres - MIN) / (MAX - MIN)
)

MAX = max(these_coords$pxl_row_in_fullres)
MIN = min(these_coords$pxl_row_in_fullres)
these_coords$array_col_new = round(
    (NUM_COL - 1) * (these_coords$pxl_row_in_fullres - MIN) / (MAX - MIN)
)
