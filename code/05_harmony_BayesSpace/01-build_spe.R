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
rownames(colData(spe)) = rownames(spatialCoords(spe)) # tibbles lose rownames

if (any(is.na(spe$array_row_transformed))) {
    stop("Some NA transformed spot coords: number of transformed spot coords likely doesn't match ncol(spe)")
}

################################################################################
#   Use transformed spatial coordinates by default
################################################################################

message("Using transformed spatialCoords by default")

#   Swap out original spot pixel coordinates for transformed ones
spe$pxl_col_in_fullres = unname(spatialCoords(spe)[, 'pxl_col_in_fullres'])
spe$pxl_row_in_fullres = unname(spatialCoords(spe)[, 'pxl_row_in_fullres'])

spatial_coords = as.matrix(
    colData(spe)[,c('pxl_col_in_fullres_transformed', 'pxl_row_in_fullres_transformed')]
)
colnames(spatial_coords) = c('pxl_col_in_fullres', 'pxl_row_in_fullres')
rownames(spatial_coords) = rownames(spatialCoords(spe))
spatialCoords(spe) = spatial_coords

#   Swap out original spot array coordinates for transformed ones
spe$array_row_original = spe$array_row
spe$array_col_original = spe$array_col
spe$array_row = spe$array_row_transformed
spe$array_col = spe$array_col_transformed

################################################################################
#   Use merged images (one image per donor) instead of single-capture-area
#   images
################################################################################

message("Overwriting imgData(spe) with merged images (one per donor)")

img_data = readImgData(
    path = file.path(transformed_dir, unique(sample_info$donor)),
    sample_id = unique(sample_info$donor),
    imageSources = file.path(
        transformed_dir, unique(sample_info$donor), "tissue_lowres_image.png"
    ),
    scaleFactors = file.path(
        transformed_dir, unique(sample_info$donor), "scalefactors_json.json"
    ),
    load = TRUE
)

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

spe = SpatialExperiment(
    assays = assays(spe),
    reducedDims = reducedDims(spe),
    rowData = rowData(spe),
    colData = colData_fixed,
    spatialCoords = spatialCoords(spe),
    scaleFactors = scaleFactors(spe),
    imgData = img_data
)

#   Drop spots with 0 counts for all genes, and drop genes with 0 counts in
#   every spot
spe <- spe[
    rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0
]

message("Saving spe")
saveRDS(spe, out_path)

session_info()
