library(visiumStitched)
library(tidyverse)
library(here)
library(rjson)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)

info_in_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
info_out_path = here(
    'processed-data', '15_samui_imagej_comparison', 'sample_info.csv'
)
coords_dir = here(
    'processed-data', '15_samui_imagej_comparison', 'spe_inputs'
)
spe_out_dir = here(
    'processed-data', '15_samui_imagej_comparison', 'spe'
)
spe_in_dir = here(
    'processed-data', '05_harmony_BayesSpace', '03-filter_normalize_spe',
    'spe_filtered_hdf5'
)

sample_info = read_csv(info_in_path, show_col_types = FALSE) |>
    filter(`In analysis`) |>
    dplyr::rename(
        capture_area = `...1`,
        group = Brain
    ) |>
    mutate(
        imagej_xml_path = here(`XML file name`),
        imagej_image_path = str_replace(
            imagej_xml_path,
            'transformations/(.*?)\\.xml$',
            'stitched_images/\\1.png'
        )
    ) |>
    select(
        capture_area, group, spaceranger_dir, imagej_xml_path, imagej_image_path
    )

#   Read in high-resolution scale factors
sample_info$tissue_hires_scalef <- sapply(
    sample_info$spaceranger_dir,
    function(x) {
        fromJSON(file = file.path(x, "scalefactors_json.json"))$tissue_hires_scalef
    }
) |>
    unname()

sample_info = sample_info |>
    #   Keep donors that have the same scalefactors in all capture areas. This
    #   avoids potentially using donors downstream affected by the inconsistent-
    #   scaling bug between ImageJ and Samui
    group_by(group) |>
    filter(all(tissue_hires_scalef == max(tissue_hires_scalef))) |>
    #   For consistency with the workflow already done in the NAc, set
    #   intra_group_scalar to 1
    dplyr::rename(group_hires_scalef = tissue_hires_scalef) |>
    mutate(intra_group_scalar = 1)

#   Write transformed image, scalefactors, and coords to a temporary directory
prep_imagej_coords(sample_info, coords_dir)
prep_imagej_image(sample_info, coords_dir)

spe_imagej = build_spe(sample_info, coords_dir = coords_dir)
colnames(spe_imagej) = spe_imagej$key

spe_filtered = loadHDF5SummarizedExperiment(spe_in_dir)

#   Subset both objects to their shared spots. This essentially means spots that
#   passed QC (from the filtered object) and are in donors with the same scale
#   factors in all capture areas (avoiding the scaling bug)
common_spots = intersect(spe_filtered$key, spe_imagej$key)
spe_filtered = spe_filtered[, common_spots]
spe_imagej = spe_imagej[, common_spots]

#   Reduce the size of the object, since we'll only need the raw counts and
#   Harmony reducedDims for BayesSpace and PRECAST
assays(spe_filtered) = list(counts = assays(spe_filtered)$counts)

#   Include both versions of the array coordinates, the only spatial information
#   PRECAST or BayesSpace uses when clustering
temp = colData(spe_filtered) |>
    as_tibble() |>
    dplyr::rename(
        array_row_samui = array_row, array_col_samui = array_col
    ) |>
    mutate(
        array_row_imagej = spe_imagej$array_row,
        array_col_imagej = spe_imagej$array_col
    ) |>
    DataFrame()
rownames(temp) = colnames(spe_filtered)
colData(spe_filtered) = temp

saveHDF5SummarizedExperiment(spe_filtered, spe_out_dir, replace = TRUE)

session_info()
