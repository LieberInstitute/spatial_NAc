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
    file.path('tissue_positions_list.csv')

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

session_info()
