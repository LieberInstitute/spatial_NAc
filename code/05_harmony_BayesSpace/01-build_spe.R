library(spatialLIBD)
library(here)
library(tidyverse)
library(jaffelab)

sample_info_path = here('raw-data', 'sample_key_spatial_NAc.csv')
sample_info_path2 = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
transformed_dirs = here('processed-data', '04_VisiumStitcher', '{}')

#   Read in the two sources of sample info and merge
sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(sample_id = Slide) |>
    select(c(sample_id, Age, Sex, Diagnosis, Refined.transforms))
sample_info = read.csv(sample_info_path2) |>
    as_tibble() |>
    rename(sample_id = X) |>
    right_join(sample_info, by = "sample_id") |>
    filter(Refined.transforms == "Yes") |>
    mutate(spaceranger_dir = dirname(normalizePath(spaceranger_dir)))

#   Build the SPE object using the individual capture areas as samples, to start
spe = read10xVisiumWrapper(
    samples = sample_info$spaceranger_dir,
    sample_id = sample_info$sample_id,
    type = "sparse",
    data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    verbose = TRUE
)

#   Merge all transformed spot coordinates into a single tibble
coords_list = list()
for (donor in unique(sample_info$Brain)) {
    #   Read in and clean up transformed spot coordinates
    coords_list[[donor]] = file.path(transformed_dirs, 'tissue_positions_list.csv') |>
        str_replace('\\{\\}', donor) |>
        read.csv() |>
        as_tibble() |>
        mutate(
            #   Reformat to be [barcode]_[sample_id] (matching spe$key)
            key = paste(
                ss(barcode, '_[ABCD]1', 2),
                barcode |> str_replace('(.*_[ABCD]1).*', '\\1'),
                sep = '_'
            )
        ) |>
        rename(
            array_row_transformed = array_row,
            array_col_transformed = array_col,
            pxl_row_in_fullres_transformed = pxl_row_in_fullres,
            pxl_col_in_fullres_transformed = pxl_col_in_fullres
        ) |>
        select(-c(in_tissue, barcode))
}
coords = do.call(rbind, coords_list)

#   Add transformed spot coordinates to colData
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(coords, by = "key") |>
    DataFrame()

if (any(is.na(spe$array_row_transformed))) {
    stop("Some NA transformed spot coords: number of transformed spot coords likely doesn't match ncol(spe)")
}
