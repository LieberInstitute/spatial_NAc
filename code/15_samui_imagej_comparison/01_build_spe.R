library(visiumStitched)
library(tidyverse)
library(here)
library(rjson)
library(SpatialExperiment)

info_in_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
info_out_path = here(
    'processed-data', '15_samui_imagej_comparison', 'sample_info.csv'
)
coords_dir = here(
    'processed-data', '15_samui_imagej_comparison', 'spe_inputs'
)
spe_out_path = here(
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
)

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

spe = build_spe(sample_info, coords_dir = coords_dir)

#   Find overlaps and determine which spots to drop for plotting based on total
#   UMI counts
spe$sum_umi = colSums(assays(spe)$counts)
spe = add_overlap_info(spe, "sum_umi")
