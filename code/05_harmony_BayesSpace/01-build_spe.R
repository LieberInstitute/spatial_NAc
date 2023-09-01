library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(scran)
library(BiocParallel)

sample_info_path = here('raw-data', 'sample_key_spatial_NAc.csv')
sample_info_path2 = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)
transformed_dir = here('processed-data', '04_VisiumStitcher')
raw_out_path = here('processed-data', '05_harmony_BayesSpace', 'spe_raw.rds')
filtered_out_path = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
num_cores = 4

################################################################################
#   Read in the two sources of sample info and merge
################################################################################

message("Gathering sample info")
sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(sample_id = Slide) |>
    select(c(sample_id, Age, Sex, Diagnosis, In.analysis, Refined.transforms))

sample_info = read.csv(sample_info_path2) |>
    as_tibble() |>
    rename(sample_id = X) |>
    right_join(sample_info, by = "sample_id") |>
    filter(In.analysis == "Yes") |>
    mutate(spaceranger_dir = dirname(normalizePath(spaceranger_dir))) |>
    rename(slide_num = Slide.., array_num = Array.., donor = Brain)

#   Grab donors that are in analysis but don't have refined transforms
unrefined_donors = sample_info |>
    filter(In.analysis == "Yes", Refined.transforms == "No") |>
    pull(donor) |>
    unique()

#   Similarly for those that do have refined transforms
refined_donors = sample_info |>
    filter(In.analysis == "Yes", Refined.transforms == "Yes") |>
    pull(donor) |>
    unique()

#   Ensure each unrefined donor consists of just one sample
n_unrefined_samples = sample_info |>
    filter(donor %in% unrefined_donors) |>
    nrow()
stopifnot(n_unrefined_samples == length(unrefined_donors))

#   We only need certain columns in the colData
sample_info = sample_info |>
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
    data = "raw",
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
for (donor in refined_donors) {
    #   Read in and clean up transformed spot coordinates
    spot_coords = file.path(
            transformed_dir, donor, 'tissue_positions_list.csv'
        ) |>
        read.csv() |>
        as_tibble() |>
        rename(
            key = barcode,
            array_row_transformed = array_row,
            array_col_transformed = array_col,
            pxl_row_in_fullres_transformed = pxl_row_in_fullres,
            pxl_col_in_fullres_transformed = pxl_col_in_fullres
        ) |>
        select(-in_tissue)
    
    #   Read in and clean Visium Stitcher info about spot overlaps
    coords_list[[donor]] = file.path(
            transformed_dir, donor, 'vs_stitcher.csv'
        ) |>
        read.csv() |>
        as_tibble() |>
        rename(key = X) |>
        select(c(key, overlap_slide, exclude_overlapping)) |>
        left_join(spot_coords, by = "key") |>
        mutate(
            #   Reformat to be [barcode]_[sample_id] (matching spe$key)
            key = paste(
                ss(key, '_[ABCD]1', 2),
                key |> str_replace('(.*_[ABCD]1).*', '\\1'),
                sep = '_'
            )
        )
}
coords = do.call(rbind, coords_list)

#   Add transformed spot coordinates, spot-overlap info, and sample_info to
#   colData
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(coords, by = "key") |>
    left_join(sample_info, by = 'sample_id') |>
    DataFrame()
rownames(colData(spe)) = rownames(spatialCoords(spe)) # tibbles lose rownames

################################################################################
#   Use transformed spatial coordinates by default
################################################################################

message("Using transformed spatialCoords by default")

#   Swap out original spot pixel coordinates for transformed ones, where they
#   exist
spe$pxl_col_in_fullres_original = unname(spatialCoords(spe)[, 'pxl_col_in_fullres'])
spe$pxl_row_in_fullres_original = unname(spatialCoords(spe)[, 'pxl_row_in_fullres'])

coords_col = ifelse(
    is.na(spe$pxl_col_in_fullres_transformed),
    spe$pxl_col_in_fullres_original, spe$pxl_col_in_fullres_transformed
)
coords_row = ifelse(
    is.na(spe$pxl_row_in_fullres_transformed),
    spe$pxl_row_in_fullres_original, spe$pxl_row_in_fullres_transformed
)
spatialCoords(spe) = matrix(
    c(coords_col, coords_row),
    ncol = 2,
    dimnames = list(
        rownames(spatialCoords(spe)),
        c('pxl_col_in_fullres', 'pxl_row_in_fullres')
    )
)

#   Swap out original spot array coordinates for transformed ones, where they
#   exist
spe$array_row_original = spe$array_row
spe$array_col_original = spe$array_col
spe$array_row = ifelse(
    is.na(spe$array_row_transformed),
    spe$array_row_original, spe$array_row_transformed
)
spe$array_col = ifelse(
    is.na(spe$array_col_transformed),
    spe$array_col_original, spe$array_col_transformed
)

################################################################################
#   Use merged images (one image per donor) instead of single-capture-area
#   images
################################################################################

message("Overwriting imgData(spe) with merged images (one per donor)")

#   Read in the merged images for donors where they exist
img_data = readImgData(
    path = file.path(transformed_dir, refined_donors),
    sample_id = refined_donors,
    imageSources = file.path(
        transformed_dir, refined_donors, "tissue_lowres_image.png"
    ),
    scaleFactors = file.path(
        transformed_dir, refined_donors, "scalefactors_json.json"
    ),
    load = TRUE
)

#   For donors with one sample, use the images already read in from spaceranger
unrefined_ids = sample_info[[
    match(unrefined_donors, sample_info$donor), 'sample_id'
]]
img_data_unrefined = imgData(spe)[imgData(spe)$sample_id %in% unrefined_ids,]
img_data_unrefined$sample_id = sample_info$donor[
    match(img_data_unrefined$sample_id, sample_info$sample_id)
]

img_data = rbind(img_data, img_data_unrefined)

################################################################################
#   Change from [slide]_[array] to donor as sample ID
################################################################################

message("Using donor as sample ID instead of [slide]_[array]")

#   Ideally, we'd include images for individual capture areas as well as by
#   donor (from the stitching process). Due to SpatialExperiment validity checks
#   that throw errors with basic manipulations of an object, this isn't possible
#   (see https://github.com/drighelli/SpatialExperiment/blob/d1934540f5f33a0aa1ae4f886f3e0ef390c210c9/R/Validity.R#L36-L40).
#   Thus, we have to pick one sample ID variable, so we'll use donor

#   There are also very robust and strict checks preventing direct modification
#   of spe$sample_id, and breaking things if imgData sample IDs don't properly
#   match sample IDs in colData. It turns out so difficult that it appears
#   necessary to construct an entirely new SpatialExperiment (simultaneously
#   defining spe$sample_id and imgData sample IDs in a compatible manner),
#   to functionally modify these parts of the object.
#  
#   See https://github.com/drighelli/SpatialExperiment/blob/d1934540f5f33a0aa1ae4f886f3e0ef390c210c9/R/SpatialExperiment-colData.R#L76-L81
#   for the colData checks
colData_fixed = colData(spe)
colData_fixed$sample_id_original = colData_fixed$sample_id
colData_fixed$sample_id = colData_fixed$donor

#   Remove the 'X' added to colnames that start with an integer when converting
#   to a tibble at various points in the script
colnames(colData_fixed) = sub('^X([0-9])', '\\1', colnames(colData_fixed))

#   Fix data types in colData: ensure categorical data are factors
categorical_cols = c(
    'sample_id', 'donor', 'slide_num', 'array_num', 'Sex', 'Diagnosis',
    'sample_id_original', 'overlap_slide',
    colnames(colData_fixed)[grep('^10x_', colnames(colData_fixed))]
)
for (categorical_col in categorical_cols) {
    colData_fixed[, categorical_col] = as.factor(
        colData_fixed[, categorical_col]
    )
}

spe = SpatialExperiment(
    assays = assays(spe),
    reducedDims = reducedDims(spe),
    rowData = rowData(spe),
    colData = colData_fixed,
    spatialCoords = spatialCoords(spe),
    scaleFactors = scaleFactors(spe),
    imgData = img_data
)

#   Make spot IDs unique
colnames(spe) = spe$key

#   Save the full object (all spots)
message("Saving raw spe")
saveRDS(spe, raw_out_path)

################################################################################
#   Compute log-normalized counts
################################################################################

#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
spe <- spe[
    rowSums(assays(spe)$counts) > 0,
    spe$in_tissue & (colSums(assays(spe)$counts) > 0)
]

message("Running quickCluster()")

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id_original,
    block.BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

message("Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(
    spe,
    clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

table(spe$scran_quick_cluster)

message("Running checking sizeFactors()")
summary(sizeFactors(spe))

message("Running logNormCounts()")
spe <- logNormCounts(spe)

message("Saving filtered spe")
saveRDS(spe, filtered_out_path)

session_info()
