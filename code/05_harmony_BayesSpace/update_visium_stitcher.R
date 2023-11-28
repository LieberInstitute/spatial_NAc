#   The SpatialExperiments originally built in 01-build_spe.* and
#   02-filter_normalize_spe.* are great except for the VisiumSticher-related
#   colData columns. Since it would be overly time-consuming to recreate the
#   full objects from scratch, this script updates the existing objects in place
#   with the new Visium Stitcher results.

library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(jaffelab)
library(tidyverse)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)

raw_path = here('processed-data', '05_harmony_BayesSpace', 'spe_raw.rds')
filtered_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
transformed_dir = here('processed-data', '04_VisiumStitcher')
plot_dir = here('plots', '05_harmony_BayesSpace')

spe_raw = readRDS(raw_path)

################################################################################
#   Functions
################################################################################

#   Add the columns of 'coords' to colData(spe), joining on 'key'. Return the
#   updated colData
add_to_coldata = function(spe, coords) {
    orig = colData(spe)

    a = colData(spe) |>
        as_tibble() |>
        select(-c(exclude_overlapping, overlap_slide)) |>
        left_join(coords, by = "key") |>
        DataFrame()
    rownames(a) = a$key

    #   Remove the 'X' added to colnames that start with an integer when converting
    #   to a tibble
    colnames(a) = sub('^X([0-9])', '\\1', colnames(a))

    #   Ensure colData has not changed except the 'exclude_overlapping' and
    #   'overlap_slide' columns
    after = a
    after$exclude_overlapping = NULL
    after$overlap_slide = NULL
    orig$exclude_overlapping = NULL
    orig$overlap_slide = NULL
    stopifnot(identical(orig, after))

    return(a)
}

################################################################################
#   Read in Visium Stitcher results
################################################################################

coords_list = list()
for (this_donor in unique(spe_raw$donor)) {
    #   Read in spots and 'exclude_overlapping' info from Visium Stitcher. In
    #   general, spots here are a subset of those from the Samui workflow
    coords_list[[this_donor]] = file.path(
            transformed_dir, this_donor, 'vs_stitcher.csv'
        ) |>
        read.csv() |>
        as_tibble() |>
        rename(key = X) |>
        #   Reformat to be [barcode]_[sample_id] (matching spe$key)
        mutate(
            key = paste(
                ss(key, '_[ABCD]1', 2),
                key |> str_replace('(.*_[ABCD]1).*', '\\1'),
                sep = '_'
            ),
            #   Change from character to logical
            exclude_overlapping = exclude_overlapping == 'True'
        ) |>
        select(c(key, overlap_slide, exclude_overlapping))
    
    #   As a sanity check, confirm spots are identical between the in-tissue
    #   SPE and Visium Stitcher
    spe_keys = colData(spe_raw) |>
        as_tibble() |>
        filter(donor == this_donor, in_tissue) |>
        pull(key)
    stopifnot(all(sort(spe_keys) == sort(coords_list[[this_donor]]$key)))
}

coords = do.call(rbind, coords_list)

################################################################################
#   Add Visium Stitcher columns to spe_raw colData. Then resave in place
################################################################################

colData(spe_raw) = add_to_coldata(spe_raw, coords)
saveRDS(spe_raw, raw_path)

################################################################################
#   Add Visium Stitcher columns to spe_filtered colData. Then resave in place
################################################################################

#-------------------------------------------------------------------------------
#   Load in object and update colData
#-------------------------------------------------------------------------------

spe_filtered = loadHDF5SummarizedExperiment(filtered_dir)
colData(spe_filtered) = add_to_coldata(spe_filtered, coords)
if (any(is.na(spe_filtered$exclude_overlapping))) {
    stop("'exclude_overlapping' column should have no NAs in filtered object")
}

#-------------------------------------------------------------------------------
#   Exploratory plots to make sure 'exclude_overlapping' looks reasonable
#-------------------------------------------------------------------------------

#   (Here we check a sample that had visually obvious problems with
#   'exclude_overlapping' before, but the below tests are not exhaustive)

#   Plot donor and color by capture area
plot_list = list()
plot_list[[1]] = spot_plot(
    spe_filtered,
    sample_id = 'Br3942',
    var_name = 'sample_id_original',
    is_discrete = TRUE,
    spatial = TRUE
)

#   Plot the original values of 'exclude_overlapping'. Since 'spot_plot' filters
#   excluded spots, we use a temporary column to take a look
spe_filtered$exclude_temp = spe_filtered$exclude_overlapping
spe_filtered$exclude_overlapping = FALSE
plot_list[[2]] = spot_plot(
    spe_filtered,
    sample_id = 'Br3942',
    var_name = 'exclude_temp',
    is_discrete = TRUE,
    spatial = TRUE,
    title = 'Br3942_exclude_overlapping'
)

pdf(file.path(plot_dir, 'update_visium_stitcher_checks.pdf'))
print(plot_list)
dev.off()

#-------------------------------------------------------------------------------
#   Resave in place
#-------------------------------------------------------------------------------

#   Switch back columns that were modified for plotting purposes only
spe_filtered$exclude_overlapping = spe_filtered$exclude_temp
spe_filtered$exclude_temp = NULL

quickResaveHDF5SummarizedExperiment(spe_filtered)

session_info()
