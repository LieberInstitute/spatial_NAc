#   Some work with VistoSeg has already been done. This script confirms that the
#   current in-analysis samples have not been processed through VistoSeg

library(here)
library(tidyverse)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
old_image_list = here('code', 'VistoSeg', 'code', 'VNS_list.txt')

sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    filter(In.analysis == "True") |>
    rename(sample_id = X)

message('Total samples: ', nrow(sample_info))

#   No relevant raw images appear to have been segmented
sample_info$raw_image_path |>
    str_replace('\\.tif$', '_nuclei.mat') |>
    file.exists() |>
    table()

old_image_paths = readLines(old_image_list)
old_sample_ids = paste(
    str_extract(old_image_paths, 'V[0-9]{2}U[0-9]{2}-[0-9]{3}'),
    str_extract(old_image_paths, '_([ABCD]1)\\.tif$', group = 1),
    sep = '_'
)

#   Also, none of the old to-be-segmented samples are in the relevant set now
table(old_sample_ids %in% sample_info$sample_id)
