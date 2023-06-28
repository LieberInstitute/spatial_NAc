library(here)
library(tidyverse)
library(jaffelab)
library(rjson)
library(cowplot)
library(grid)
library(gridExtra)
library(getopt)
library(sessioninfo)

spec = matrix(
    c(
        "slide", "s", 1, "character", "Slide number",
        "arrays", "a", 1, "character", "Array number"
    ),
    byrow = TRUE, ncol = 5
)
opt = getopt(spec)

#   Column names for raw tissue_positions_list.csv files
tissue_colnames = c(
    "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres",
    "pxl_col_in_fullres"
)

tissue_path = here(
    'processed-data', '02_image_stitching',
    sprintf('tissue_positions_%s_%s.csv', opt$slide, opt$arrays)
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

coords_path_out = here(
    'processed-data', '02_image_stitching',
    sprintf('spatial_coords_%s_%s.csv', opt$slide, opt$arrays)
)

plot_dir = here(
    'plots', '02_image_stitching', paste(opt$slide, opt$arrays, sep = '_')
)

dir.create(plot_dir, showWarnings = FALSE)

#   55-micrometer diameter for Visium spot; 100 micrometers between spots;
#   65-micrometer spot diameter used in 'spot_diameter_fullres' calculation for
#   spaceranger JSON. See documentation for respective quantities. The
#   difference between 55 and 65 does indeed exist and is properly documented,
#   but is likely a bug in the sense the choice was probably unintentional
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487812-What-is-the-size-of-the-spots-on-the-Visium-Gene-Expression-Slide-
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487892-How-much-space-is-there-between-spots-referred-to-as-white-space-
#   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
SPOT_DIAMETER_M = 55e-6
SPOT_DIAMETER_JSON_M = 65e-6
INTER_SPOT_DIST_M = 100e-6

set.seed(0)

################################################################################
#   Functions
################################################################################

build_array = function(max_row, max_col) {
    new_pairs = c(
        outer(
            seq(0, 2 * (max_row %/% 2) + 2, 2),
            seq(0, 2 * (max_col %/% 2) + 2, 2),
            paste
        ),
        outer(
            seq(1, 2 * (max_row %/% 2) + 1, 2),
            seq(1, 2 * (max_col %/% 2) + 1, 2),
            paste
        )
    )

    new_array = tibble(
            array_row = as.integer(ss(new_pairs, ' ', 1)),
            array_col = as.integer(ss(new_pairs, ' ', 2))
        ) |>
        #   The only exception to the array pattern is that (0, 0) cannot exist
        filter(!((array_row == 0) & (array_col == 0)))

    return(new_array)
}

#   Stop if the tibble-like 'coords', expected to contain columns 'array_row'
#   and 'array_col', represents an invalid Visium array
validate_array = function(coords) {
    #   Even array rows can only use even column indices
    coords |>
        filter(array_row %% 2 == 0) |>
        summarize(a = all(array_col %% 2 == 0)) |>
        pull(a) |>
        stopifnot()

    #   Odd array rows can only use odd column indices
    coords |>
        filter(array_row %% 2 == 1) |>
        summarize(a = all(array_col %% 2 == 1)) |>
        pull(a) |>
        stopifnot()

    #   Check lower bound of array row and col (note we're allowing arbitrary
    #   maximum values rather than the convention of 78 rows and 128 columns)
    stopifnot(min(coords$array_row) == 0)
    stopifnot(min(coords$array_col) == 0)

    #   Check an eccentric detail of Visium arrays: (0, 0) cannot exist
    if (any((coords$array_row == 0) & (coords$array_col == 0))) {
        stop("The invalid array coordinate (0, 0) exists after fitting")
    }
}

refine_fit = function(x, y, INTERVAL_X, INTERVAL_Y) {
    #   Round x to the nearest integer, and track the error from doing so in the
    #   variable 'dx'
    dx = x - clean_round(x)
    x = x - dx

    #   Given x, round y to the nearest valid integer (y must be even iff x is),
    #   and track the error from doing so in the variable 'dy'
    dy = rep(0, length(y))
    dy[x %% 2 == 0] = y[x %% 2 == 0] - clean_round(y[x %% 2 == 0] / 2) * 2
    dy[x %% 2 == 1] = y[x %% 2 == 1] - (clean_round(y[x %% 2 == 1] / 2 - 0.5) * 2 + 1)
    y = y - dy

    #   Summarize error in Euclidean distance
    error = sqrt((INTERVAL_X * dx) ** 2 + (INTERVAL_Y * dy) ** 2)
    return(list(x, y, error))
}

#   Given 'coords', a tibble whose rows represent samples of the same
#   donor, and containing columns 'array_row', 'array_col',
#   'pxl_row_in_fullres', and 'pxl_col_in_fullres', modify the 'array_row' and
#   'array_col' columns to represent a larger Visium capture area containing
#   all samples in a common coordinate system. 'inter_spot_dist_px' gives the
#   distance between any 2 spots in the new coordinates, and so the number of
#   array rows/cols generally changes from the Visium standards of 78 and 128
#   (and even changes in ratio between num rows and num cols). Return the same
#   tibble with modified 'array_row' + 'array_col' columns, as well as new
#   'pxl_row_rounded' and 'pxl_col_rounded' columns representing the pixel
#   coordinates rounded to the nearest exact array coordinates.
fit_to_array = function(coords, inter_spot_dist_px) {
    MIN_ROW = min(coords$pxl_col_in_fullres)
    INTERVAL_ROW = inter_spot_dist_px * cos(pi / 6)
        
    MIN_COL = min(coords$pxl_row_in_fullres)
    INTERVAL_COL = inter_spot_dist_px / 2

    #   Calculate what "ideal" array rows and cols should be (allowed to be any
    #   float). Don't round yet
    array_row_temp = (coords$pxl_col_in_fullres - MIN_ROW) / 
        INTERVAL_ROW

    array_col_temp = (coords$pxl_row_in_fullres - MIN_COL) /
        INTERVAL_COL

    #   For now, find the nearest row first, then round to the nearest possible
    #   column given the row
    temp = refine_fit(array_row_temp, array_col_temp, INTERVAL_ROW, INTERVAL_COL)
    error_row_first = temp[[3]]
    coords$array_row = temp[[1]]
    coords$array_col = temp[[2]]

    #   Perform the opposite order (column then row). When this ordering results
    #   in lower error, use it instead
    temp = refine_fit(array_col_temp, array_row_temp, INTERVAL_COL, INTERVAL_ROW)
    error_col_first = temp[[3]]
    coords$array_row[error_row_first > error_col_first] = temp[[2]][
        error_row_first > error_col_first
    ]
    coords$array_col[error_row_first > error_col_first] = temp[[1]][
        error_row_first > error_col_first
    ]

    #   Now make new pixel columns based on just the array values (these columns
    #   give the coordinates for given array row/cols)
    coords$pxl_col_rounded = MIN_ROW + coords$array_row * INTERVAL_ROW
    coords$pxl_row_rounded = MIN_COL + coords$array_col * INTERVAL_COL

    #-------------------------------------------------------------------------------
    #   array (0, 0) does not exist on an ordinary Visium array. Move any such
    #   values to the nearest alternatives
    #-------------------------------------------------------------------------------

    #   Nearest points to (0, 0) are (0, 2) and (1, 1):
    array_02 = c(MIN_ROW, MIN_COL + 2 * INTERVAL_COL)
    array_11 = c(MIN_ROW + INTERVAL_ROW, MIN_COL + INTERVAL_COL)

    #   Determine the distances to those nearest points
    dist_coords = coords |>
        filter(array_row == 0, array_col == 0) |>
        mutate(
            dist_02 = (pxl_col_in_fullres - array_02[1]) ** 2 +
                (pxl_row_in_fullres - array_02[2]) ** 2,
            dist_11 = (pxl_col_in_fullres - array_11[1]) ** 2 +
                (pxl_row_in_fullres - array_11[2]) ** 2,
        )

    #   Move any instances of (0, 0) to the nearest alternative
    indices = (coords$array_row == 0) & (coords$array_col == 0)
    coords[indices, 'array_row'] = ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 0, 1
    )
    coords[indices, 'array_col'] = ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 2, 1
    )

    #   Verify the newly assigned array row and cols have reasonable values
    validate_array(coords)

    return(coords)
}

#   Given 'coords', a tibble whose rows represent samples of the same
#   donor, and containing columns 'array_row', 'array_col',
#   'pxl_row_in_fullres', and 'pxl_col_in_fullres', modify the 'array_row' and
#   'array_col' columns to represent a larger Visium capture area containing
#   all samples in a common coordinate system. 'inter_spot_dist_px' gives the
#   distance between any 2 spots in the new coordinates, and so the number of
#   array rows/cols generally changes from the Visium standards of 78 and 128
#   (and even changes in ratio between num rows and num cols). Return the same
#   tibble with modified 'array_row' + 'array_col' columns, as well as new
#   'pxl_row_rounded' and 'pxl_col_rounded' columns representing the pixel
#   coordinates rounded to the nearest exact array coordinates.
fit_to_array2 = function(coords, inter_spot_dist_px) {
    MIN_ROW = min(coords$pxl_col_in_fullres)
    INTERVAL_ROW = inter_spot_dist_px * cos(pi / 6)
        
    MIN_COL = min(coords$pxl_row_in_fullres)
    INTERVAL_COL = inter_spot_dist_px / 2

    #   Calculate what "ideal" array rows and cols should be (allowed to be any
    #   integer)
    coords$array_row = clean_round(
        (coords$pxl_col_in_fullres - MIN_ROW) / INTERVAL_ROW
    )

    coords$array_col = clean_round(
        (coords$pxl_row_in_fullres - MIN_COL) / INTERVAL_COL
    )

    #   Form a "new Visium array" spanning coordinates of all samples
    new_array = build_array(max(coords$array_row), max(coords$array_col))
    new_array$pxl_col_rounded = MIN_ROW + new_array$array_row * INTERVAL_ROW
    new_array$pxl_row_rounded = MIN_COL + new_array$array_col * INTERVAL_COL

    for (i in 1:nrow(coords)) {
        neighbors = new_array |>
            #   Grab the neighbors and the spot itself
            filter(
                abs(array_row - coords$array_row[i]) <= 1,
                abs(array_col - coords$array_col[i]) <= 2
            ) |>
            #   Compute the exact distance from the initial spot to its
            #   potential neighbors (actually, squared distance, to save an
            #   extra calculation)
            mutate(
                dist = (pxl_row_rounded - coords$pxl_row_in_fullres[i]) ** 2 +
                    (pxl_col_rounded - coords$pxl_col_in_fullres[i]) ** 2
            )
        
        #   Now overwrite the approximate array-coordinate estimate with the
        #   precise one that uses Euclidean distance (in pixels)
        coords[i, c('array_row', 'array_col')] = neighbors[
            which.min(neighbors$dist), c('array_row', 'array_col')
        ]
    }

    #   Now make new pixel columns based on just the array values (these columns
    #   give the coordinates for given array row/cols)
    coords$pxl_col_rounded = MIN_ROW + coords$array_row * INTERVAL_ROW
    coords$pxl_row_rounded = MIN_COL + coords$array_col * INTERVAL_COL

    #-------------------------------------------------------------------------------
    #   array (0, 0) does not exist on an ordinary Visium array. Move any such
    #   values to the nearest alternatives
    #-------------------------------------------------------------------------------

    #   Nearest points to (0, 0) are (0, 2) and (1, 1):
    array_02 = c(MIN_ROW, MIN_COL + 2 * INTERVAL_COL)
    array_11 = c(MIN_ROW + INTERVAL_ROW, MIN_COL + INTERVAL_COL)

    #   Determine the distances to those nearest points
    dist_coords = coords |>
        filter(array_row == 0, array_col == 0) |>
        mutate(
            dist_02 = (pxl_col_in_fullres - array_02[1]) ** 2 +
                (pxl_row_in_fullres - array_02[2]) ** 2,
            dist_11 = (pxl_col_in_fullres - array_11[1]) ** 2 +
                (pxl_row_in_fullres - array_11[2]) ** 2,
        )

    #   Move any instances of (0, 0) to the nearest alternative
    indices = (coords$array_row == 0) & (coords$array_col == 0)
    coords[indices, 'array_row'] = ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 0, 1
    )
    coords[indices, 'array_col'] = ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 2, 1
    )

    #   Verify the newly assigned array row and cols have reasonable values
    validate_array(coords)

    return(coords)
}

#   Round to the nearest integer, always rounding up at 0.5. This consistent
#   behavior is favorable for our application, where we want to minimize
#   duplicate mappings
clean_round = function(x) {
    return(floor(x) + ((x * 10) %% 10 >= 5))
}

################################################################################
#   Read in sample info and spot coordinate info
################################################################################

coords = read.csv(tissue_path) |>
    as_tibble() |>
    mutate(sample_id = paste(ss(key, '_', 2), ss(key, '_', 3), sep = '_'))

sample_ids = paste(opt$slide, strsplit(opt$arrays, '_')[[1]], sep = '_')
sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(
        sample_id = X,
        array_num = 'Array..',
        slide_num = 'Slide..' 
    ) |>
    filter(sample_id %in% sample_ids)

#-------------------------------------------------------------------------------
#   Read in the untransformed and unfiltered (by "in tissue") spot coordinates
#   to verify the relationship between array and pixel coordinates
#-------------------------------------------------------------------------------

raw_tissue_paths = sample_info |>
    pull(spaceranger_dir) |>
    file.path('tissue_positions_list.csv')

#   Read in coordinates, keeping track of sample ID as well
raw_tissue_list = list()
for (raw_tissue_path in raw_tissue_paths) {
    raw_tissue_list[[tissue_path]] = read.csv(raw_tissue_path, col.names = tissue_colnames)
    raw_tissue_list[[tissue_path]]$sample_id = tissue_path |>
        str_extract('/[^/]*_[ABCD]1/') |>
        str_replace_all('/', '')
}
raw_coords = do.call(rbind, raw_tissue_list) |> as_tibble()

print("Correlation relationship between array (a) and pixel (p) row/col:")
cor_raw_coords = raw_coords |>
    summarize(
        arow_v_prow = cor(array_row, pxl_row_in_fullres),
        arow_v_pcol = cor(array_row, pxl_col_in_fullres),
        acol_v_prow = cor(array_col, pxl_row_in_fullres),
        acol_v_pcol = cor(array_col, pxl_col_in_fullres)
    )
print(cor_raw_coords)

#   We expect array_row to "line up with" pixel_col
tol = 0.001
stopifnot(all(abs(cor_raw_coords$arow_v_prow) < tol))
stopifnot(all(1 - abs(cor_raw_coords$arow_v_pcol) < tol))

#-------------------------------------------------------------------------------
#   Select an appropriate number of array rows and columns for this slide such
#   that the distance between spots on the new grid matches the Visium standard
#-------------------------------------------------------------------------------

sr_json = fromJSON(
    file = file.path(
        sample_info$spaceranger_dir[1], "scalefactors_json.json"
    )
)

#   Since we have the spot diameter both in pixels and meters, compute the
#   image's pixel/m ratio. Then use that to compute the distance between spots
#   in pixels
PX_PER_M = sr_json$spot_diameter_fullres / SPOT_DIAMETER_JSON_M
INTER_SPOT_DIST_PX = INTER_SPOT_DIST_M * PX_PER_M

#   Sanity check on existing coordinates: spots should be 100um apart. If not,
#   modify. TODO: look for pairs of array coordinates that exist; don't hardcode
#   10, 11
a = coords |>
    filter(array_row == 10, array_col == 80) |>
    head(1) |>
    select(pxl_row_in_fullres, pxl_col_in_fullres)

b = coords |>
    filter(array_row == 11, array_col == 81) |>
    head(1) |>
    select(pxl_row_in_fullres, pxl_col_in_fullres)

tol = 2
observed_dist = sqrt(rowSums((a - b) ** 2))
if (abs(observed_dist - INTER_SPOT_DIST_PX) > tol) {
    stop('Transformed spots output from ImageJ/Samui were not close to 100um apart!')
}

################################################################################
#   Visually assess array_row and array_col after fitting to "new grid"
################################################################################

#   Adjust 'array_row' and 'array_col' with values appropriate for the new
#   coordinate system (a larger Visium grid with equal inter-spot distances)
coords$pxl_row_in_fullres[coords$sample_id == "V12D07-333_A1"] = coords$pxl_row_in_fullres[coords$sample_id == "V12D07-333_A1"] + 11 * INTER_SPOT_DIST_PX / 2
# coords$pxl_col_in_fullres[coords$sample_id == "V12D07-333_A1"] = coords$pxl_col_in_fullres[coords$sample_id == "V12D07-333_A1"] + INTER_SPOT_DIST_PX / 4
coords = fit_to_array2(coords, INTER_SPOT_DIST_PX)

write.csv(
    coords |> select(!c(pxl_row_rounded, pxl_col_rounded)),
    coords_path_out,
    row.names = FALSE, quote = FALSE
)

#-------------------------------------------------------------------------------
#   Spot plots
#-------------------------------------------------------------------------------

#   Plot the transformed spots as-is for reference
p = ggplot(coords) +
    geom_point(
        aes(
            x = pxl_col_in_fullres,
            y = max(pxl_row_in_fullres) - pxl_row_in_fullres,
            color = sample_id
        ),
        size = 0.1
    ) +
    coord_fixed() +
    labs(
        x = 'Full-res pixels (width)', y = 'Full-res pixels (height)',
        color = 'Sample ID'
    ) + 
    guides(color = guide_legend(override.aes = list(size = 1.5)))

pdf(file.path(plot_dir, 'raw_spots.pdf'))
print(p)
dev.off()

#   Plot the spots aligned to the new Visium grid
p = ggplot(coords) +
    geom_point(
        aes(
            x = pxl_col_rounded,
            y = max(pxl_row_rounded) - pxl_row_rounded,
            color = sample_id
        ),
        size = 0.1
    ) +
    coord_fixed() +
    labs(
        x = 'Full-res pixels (width)', y = 'Full-res pixels (height)',
        color = 'Sample ID', title = 'Spot Positions (Aligned to Array)'
    ) + 
    guides(color = guide_legend(override.aes = list(size = 1.5)))

pdf(file.path(plot_dir, 'aligned_spots.pdf'))
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   Plots measuring mapping of multiple spots to the same array coordinates
#-------------------------------------------------------------------------------

#   Verify that for a single sample, mappings were unique
dup_coords = coords |>
    group_by(sample_id, array_col, array_row) |>
    summarize(n = n()) |>
    ungroup()

p1 = ggplot(dup_coords) +
    geom_histogram(
        aes(x = n, fill = sample_id), bins = length(unique(dup_coords$n)),
        color = "black"
    ) +
    facet_wrap(~sample_id, nrow = 1) +
    labs(
        x = "Num spots per array coordinate", fill = "Sample ID",
        title = "Within-sample duplication of spot mappings"
    )

#   Then check how many spots mapped to the same array coordinates, agnostic
#   to sample
dup_coords = coords |>
    group_by(array_col, array_row) |>
    summarize(n = n()) |>
    ungroup()

p2 = ggplot(dup_coords) +
    geom_histogram(aes(x = n), bins = length(unique(dup_coords$n))) +
    labs(
        x = "Num spots per array coordinate",
        title = "Overall duplication of spot mappings"
    )

pdf(file.path(plot_dir, 'spot_duplication.pdf'))
plot_grid(p1, p2, nrow = 2)
dev.off()

p = coords |>
    #   Clean tibble into a nice format for ggplot
    pivot_longer(
        cols = c(
            pxl_row_in_fullres, pxl_row_rounded, pxl_col_in_fullres,
            pxl_col_rounded
        ),
        values_to = "pxl_coord",
        names_pattern = "^(pxl_col|pxl_row)_(in_fullres|rounded)$",
        names_to = c("axis", "position")
    ) |>
    pivot_wider(names_from = axis, values_from = pxl_coord) |>
    mutate(
        #   Use better names
        position = case_when(
            position == "in_fullres" ~ "original",
            position == "rounded" ~ "aligned",
            TRUE ~ position
        )
    ) |>
    #   Plot the original spot coordiates on top of their alignments,
    #   faceted by sample
    ggplot() +
        geom_point(
            aes(x = pxl_col, y = max(pxl_row) - pxl_row, color = position),
            size = 0.1
        ) +
        facet_wrap(~sample_id) +
        coord_fixed()

pdf(file.path(plot_dir, 'alignment_by_sample.pdf'), width = 14, height = 14)
print(p)
dev.off()

################################################################################
#   Measure error in aligning (rounding) pixel coordinates to fit array
################################################################################

#   Mean alignment error in spot diameters
print("Mean alignment error in spot diameters:")
coords |>
    summarize(
        d = mean(
            (pxl_row_in_fullres - pxl_row_rounded) ** 2 +
            (pxl_col_in_fullres - pxl_col_rounded) ** 2
        ) ** 0.5 / sr_json$spot_diameter_fullres
    ) |>
    pull(d) |>
    round(2) |>
    print()

#-------------------------------------------------------------------------------
#   Explore why/ where duplication occurs
#-------------------------------------------------------------------------------

#   Take 4 random duplicated mappings from each sample ID
dup_df_orig = coords |>
    group_by(array_row, array_col, sample_id) |>
    mutate(num_members = n()) |>
    filter(num_members > 1) |>
    slice_head(n = 1) |>
    group_by(sample_id) |>
    sample_n(size = 4, replace = TRUE) |>
    select(!num_members)

plot_list = list()
for (i in 1:nrow(dup_df_orig)) {
    #   Grab the set of source spots mapping to the same destination spot,
    #   mark them as duplicates, and use their original pixel coordinates
    dup_df = coords |>
        filter(
            sample_id == dup_df_orig$sample_id[i],
            array_row == dup_df_orig$array_row[i],
            array_col == dup_df_orig$array_col[i]
        ) |>
        mutate(
            coord_type = "duplicate",
            pxl_row = pxl_row_in_fullres,
            pxl_col = pxl_col_in_fullres
        )
    
    #   Grab the set of nearby source spots surrounding the duplicate mapping,
    #   and use the rounded pixel coordinates from them
    dup_df = coords |>
        filter(
            sample_id == dup_df_orig$sample_id[i],
            abs(array_row - dup_df_orig$array_row[i]) <= 4,
            abs(array_col - dup_df_orig$array_col[i]) <= 7,
        ) |>
        mutate(
            coord_type = "nearby",
            pxl_row = pxl_row_rounded,
            pxl_col = pxl_col_rounded
        ) |>
        rbind(dup_df)
    
    plot_list[[i]] = ggplot(dup_df) +
        geom_point(aes(x = pxl_col, y = pxl_row, color = coord_type)) +
        coord_fixed()

    #   Grab the legend from the first plot (which matches all other legends)
    if (i == 1) {
        legend = get_legend(plot_list[[i]])
    }
    
    plot_list[[i]] = plot_list[[i]] +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none"
        )
}

plot_list[['legend']] = legend
pdf(file.path(plot_dir, 'duplicated_sites.pdf'))
grid.arrange(
    grobs = plot_list, bottom = "pixel col", left = "pixel row", ncol = 4
)
dev.off()

#   Check if spots actually mapped to the closest new spot
a = dup_df |>
    filter(coord_type == "duplicate") 

b = dup_df |>
    filter(coord_type == "nearby")

#   These shouldn't be different, but are for first sample
b[which.min(rowSums((b |> select(pxl_row_rounded, pxl_col_rounded) - as.matrix(a[2,] |> select(pxl_row_in_fullres, pxl_col_in_fullres))) ** 2)), c('array_row', 'array_col')]
a[1, c('array_row', 'array_col')]

session_info()
