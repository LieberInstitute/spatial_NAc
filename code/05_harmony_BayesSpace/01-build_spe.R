library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

sample_info_path = here('raw-data', 'sample_key_spatial_NAc.csv')
sample_info_path2 = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
transformed_dir = here('processed-data', '04_VisiumStitcher')
out_path = here('processed-data', '05_harmony_BayesSpace', 'spe.rds')

#   Read in the two sources of sample info and merge
message("Gathering sample info")
sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(sample_id = Slide) |>
    select(c(sample_id, Age, Sex, Diagnosis, Refined.transforms))

sample_info = read.csv(sample_info_path2) |>
    as_tibble() |>
    rename(sample_id = X) |>
    right_join(sample_info, by = "sample_id") |>
    filter(Refined.transforms == "Yes") |>
    mutate(spaceranger_dir = dirname(normalizePath(spaceranger_dir))) |>
    rename(slide_num = Slide.., array_num = Array.., donor = Brain) |>
    select(
        c(
            sample_id, donor, slide_num, array_num, spaceranger_dir,
            raw_image_path, Age, Sex, Diagnosis
        )
    )

################################################################################
#   Build the SPE object using [slide]_[capture area] as samples, to start
################################################################################

message("Building SpatialExperiment using [slide]_[capture area] as sample ID")
spe = read10xVisiumWrapper(
    samples = sample_info$spaceranger_dir,
    sample_id = sample_info$sample_id,
    type = "sparse",
    data = "filtered",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    verbose = TRUE
)

################################################################################
#   Read in transformed spot coordinates and add to colData
################################################################################

message("Adding transformed spot coordinates and sample info to colData")

#   Merge all transformed spot coordinates into a single tibble
coords_list = list()
for (donor in unique(sample_info$donor)) {
    #   Read in and clean up transformed spot coordinates
    coords_list[[donor]] = file.path(
            transformed_dir, donor, 'tissue_positions_list.csv'
        ) |>
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

#   Add transformed spot coordinates to colData. Also add sample_info
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(coords, by = "key") |>
    left_join(sample_info, by = 'sample_id') |>
    DataFrame()

if (any(is.na(spe$array_row_transformed))) {
    stop("Some NA transformed spot coords: number of transformed spot coords likely doesn't match ncol(spe)")
}

################################################################################
#   Add merged images (one image per donor)
################################################################################

message("Adding merged images (one per donor) to imgData(spe)")

merged_images = readImgData(
    path = file.path(transformed_dir, unique(sample_info$donor)),
    sample_id = unique(sample_info$donor),
    imageSources = file.path(
        transformed_dir, unique(sample_info$donor), "tissue_hires_image.png"
    ),
    scaleFactors = file.path(
        transformed_dir, unique(sample_info$donor), "scalefactors_json.json"
    ),
    load = TRUE
)

imgData(spe) = rbind(imgData(spe), merged_images)

message("Saving spe")
saveRDS(spe, out_path)

session_info()
