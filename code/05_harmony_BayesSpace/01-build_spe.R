library(spatialLIBD)
library(here)
library(tidyverse)

sample_info_path = here('raw-data', 'sample_key_spatial_NAc.csv')
sample_info_path2 = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

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
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    verbose = TRUE
)
