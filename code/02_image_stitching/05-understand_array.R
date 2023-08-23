library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

tissue_colnames = c(
    "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres",
    "pxl_col_in_fullres"
)

sample_info_path = here(
    'processed-data', '02_image_stitching', 'sample_info_clean.csv'
)

plot_path = here(
    'plots', '02_image_stitching', 'array_vs_pixel.pdf'
)

dir.create(dirname(plot_path), showWarnings = FALSE)

#   Read in sample_info and use better colnames
sample_info = read.csv(sample_info_path) |>
    as_tibble() |>
    rename(
        sample_id = X,
        array_num = 'Array..',
        slide_num = 'Slide..' 
    )

#   Get the paths to the original tissue positions for several example samples
#   in this dataset
tissue_paths = sample_info |>
    filter(
        slide_num %in% c("V11U08-082", "V11U23-406", "V12D07-074"),
        array_num == 'A1'
    ) |>
    pull(spaceranger_dir) |>
    list.files(pattern = '^tissue_positions(_list|)\\.csv$', full.names = TRUE)

#   Read in coordinates, keeping track of sample ID as well
tissue_list = list()
for (tissue_path in tissue_paths) {
    tissue_list[[tissue_path]] = read.csv(tissue_path, col.names = tissue_colnames)
    tissue_list[[tissue_path]]$sample_id = tissue_path |>
        str_extract('/[^/]*_A1/') |>
        str_replace_all('/', '')
}
coords = do.call(rbind, tissue_list) |> as_tibble()

#   Plot the relationship between pixel and array coordinates, which seems to
#   vary widely between samples
p = ggplot(coords) +
    geom_point(
        aes(
            x = pxl_row_in_fullres, y = pxl_col_in_fullres,
            color = array_row %% 10,
            alpha = 0.6 + 0.4 * array_col / max(array_col)
        )
    ) +
    coord_fixed() +
    facet_wrap(~sample_id) +
    scale_color_continuous(type = 'viridis')

pdf(plot_path, width = 10, height = 5)
print(p)
dev.off()

#   Given spatial coords 'coords', get the distance between neighboring spots
#   to the spot at array value (start_row, start_col)
get_dist = function(start_row, start_col, coords) {
    a = coords |>
        group_by(sample_id) |>
        filter(array_row == start_row, array_col == start_col) |>
        ungroup() |>
        select(pxl_row_in_fullres, pxl_col_in_fullres)
    b = coords |>
        group_by(sample_id) |>
        filter(array_row == start_row, array_col == start_col + 2) |>
        ungroup() |>
        select(pxl_row_in_fullres, pxl_col_in_fullres)

    observed_dist = sqrt(rowSums((a - b) ** 2))
    return(observed_dist)
}

#   Try 25 random pairs of coordinates across all samples, to verify they are
#   consistent within and between arrays. Figured out the slight variation is
#   due to rounding spot centers to the nearest integer pixel in full resolution
num_points = 25
rows_vec = sample(1:38, num_points, replace = FALSE) * 2
cols_vec = sample(0:62, num_points, replace = FALSE) * 2
d = map2(rows_vec, cols_vec, get_dist, coords) |>
    unlist() |>
    matrix(
        ncol = coords |> pull(sample_id) |> unique() |> length()
    )
print(d)

################################################################################
#   Simulate data to verify if strange artifacting when aligning spots to a
#   "new Visium grid" still occurs
################################################################################

build_array = function(max_row, max_col, spot_distance = 1) {
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
        filter(!((array_row == 0) & (array_col == 0))) |>
        mutate(
            pxl_row = array_col * spot_distance / 2,
            pxl_col = array_row * spot_distance * cos(pi / 6)
        )       

    return(new_array)
}

rotate = function(coords, theta) {
    rot = matrix(c(cos(theta), sin(theta), -1 * sin(theta), cos(theta)), nrow = 2)
    return(t(rot %*% t(as.matrix(coords[, c('pxl_row', 'pxl_col')]))))
}

new_array = build_array(20, 40) |> mutate(coord_type = "new")
old_array = build_array(20, 40) |>
    mutate(
        pxl_row = 0.99 * pxl_row + 0.49,
        pxl_col = 0.995 * pxl_col,
        coord_type = "old"
    )
theta = 0 # 2 * pi / 180
old_array[, c('pxl_row', 'pxl_col')] = rotate(old_array[, c('pxl_row', 'pxl_col')], theta)

#   Verify we actually have valid-looking Visium arrays
ggplot(rbind(new_array, old_array)) +
    geom_point(
        aes(x = pxl_col, y = max(pxl_row) - pxl_row, color = coord_type)
    ) +
    coord_fixed()

old_array$pxl_row_rounded = NA
old_array$pxl_col_rounded = NA
for (i in 1:nrow(old_array)) {
    dist_sq = new_array |>
        mutate(
            dist = 0.5 * (pxl_row - old_array$pxl_row[i]) ** 2 +
                cos(pi / 6) * (pxl_col - old_array$pxl_col[i]) ** 2
        ) |>
        pull(dist)
    
    old_array[i, c('pxl_row_rounded', 'pxl_col_rounded')] = new_array[
        which.min(dist_sq), c('pxl_row', 'pxl_col')
    ]
}

#   Artifacting does not seem to occur with any combination of slight scaling
#   and rotating
old_array |>
    #   Clean tibble into a nice format for ggplot
    pivot_longer(
        cols = c(
            pxl_row, pxl_row_rounded, pxl_col, pxl_col_rounded
        ),
        values_to = "pxl_coord",
        names_pattern = "^(pxl_col|pxl_row)(.*)",
        names_to = c("axis", "position")
    ) |>
    pivot_wider(names_from = axis, values_from = pxl_coord) |>
    mutate(
        #   Use better names
        position = case_when(
            position == "" ~ "original",
            position == "rounded" ~ "aligned",
            TRUE ~ position
        )
    ) |>
    filter(position == "_rounded") |>
    #   Plot the original spot coordiates on top of their alignments
    ggplot() +
        geom_point(
            aes(x = pxl_col, y = max(pxl_row) - pxl_row, color = position),
        ) +
        coord_fixed()

#   Investigated the spaceranger JSONs. Spot diameters are not exactly the same
#   number of pixels across samples, but they are within 1%
json_paths = sample_info |>
    pull(spaceranger_dir) |>
    file.path('scalefactors_json.json')

json_list = list()
for (json_path in json_paths) {
    json_list[[json_path]] = fromJSON(file = json_path)
}

#   This was the problematic sample (with artifacting)
this_coords = coords |>
    filter(sample_id == "V12D07-333_A1")

#   Given a vector of array_row and array_col, return a logical vector that is
#   TRUE where those array coordinates are not present in 'this_coords'
empty_fun = function(row_vec, col_vec) {
    is_empty = c()
    for (i in 1:length(row_vec)) {
        num_matches = this_coords |>
            filter(array_row == row_vec[i], array_col == col_vec[i]) |>
            nrow()
        is_empty = c(is_empty, num_matches == 0)
    }

    return(is_empty)
}

#   Using 'new_array' defined in 06-spatial_coords.R, form a tibble of array
#   coordinates that have no spots mapping to them
missed_array = new_array |>
    #   Just the array overlapping this_coords
    filter(
        array_row >= min(this_coords$array_row),
        array_row <= max(this_coords$array_row),
        array_col >= min(this_coords$array_col),
        array_col <= max(this_coords$array_col),
        empty_fun(array_row, array_col)
    )

#   No matter which row of 'missing_array' is used, the mapping appears to be
#   correct (minimize distance)
this_coords |>
    filter(
        abs(array_row - missed_array$array_row[101]) <= 1,
        abs(array_col - missed_array$array_col[101]) <= 2
    ) |>
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
    ggplot() +
        geom_point(
            aes(
                x = pxl_col,
                y = max(pxl_row) - pxl_row,
                color = position
            )
        ) +
        coord_fixed()
    
#   Verify we're grabbing the correct spots (those with no mappings to them)
new_array |>
    #   Just the array overlapping this_coords
    filter(
        array_row >= min(this_coords$array_row),
        array_row <= max(this_coords$array_row),
        array_col >= min(this_coords$array_col),
        array_col <= max(this_coords$array_col)
    ) |>
    mutate(
        is_empty = empty_fun(array_row, array_col)
    ) |>
    ggplot() +
        geom_point(
            aes(
                x = pxl_col_rounded,
                y = max(pxl_row_rounded) - pxl_row_rounded,
                color = is_empty
            )
        ) +
        coord_fixed()

a = raw_coords |>
    filter(sample_id == "V12D07-333_A1")

#   Get the pixel distance between adjacent spots
get_spot_distance = function(index, a) {
    tibble_a = a |>
        filter(array_row == index, array_col == index) |>
        select(pxl_row_in_fullres, pxl_col_in_fullres)
    tibble_b = a |>
        filter(array_row == index+1, array_col == index+1) |>
        select(pxl_row_in_fullres, pxl_col_in_fullres)

    #   Euclidean distance in pixels
    sqrt(rowSums((tibble_b - tibble_a) ** 2))
}

#   Pixel distances do vary slightly! Must be because they must be integers by
#   definition, but attempt to mark a continuous quantity
sapply(1:20, get_spot_distance, a)
