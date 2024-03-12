library(here)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(Polychrome)
data(palette36)

precast_path <- here("processed-data", "10_precast", "PRECAST_k%s.csv")
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered_hdf5"
)
plot_dir <- here("plots", "10_precast", "figures")

#   Define color values for specific values of k
k_colors = list()
k_colors[[2]] =  palette36[c(8, 6)]
names(k_colors[[2]]) = 1:2
k_colors[[4]] = c(palette36[c(7, 8, 6)], "#45200D")
names(k_colors[[4]]) = 1:4

good_donors = c('Br8492', 'Br6522', 'Br6423')

dir.create(plot_dir, showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)

################################################################################
#   Import PRECAST results and append to colData
################################################################################

result_list <- list()
for (K in 2:28) {
    result_list[[K]] <- sprintf(precast_path, K) |>
        read_csv(show_col_types = FALSE) |>
        select(c(key, cluster)) |>
        mutate(k = K)
}
precast_results <- do.call(rbind, result_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "precast_k"
    )

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

################################################################################
#   Figure plots
################################################################################

#   Loop through values of k where the color palette is defined; plot spot
#   plots for select donors in individual PDF pages
for (k in which(!sapply(k_colors, is.null))) {
    plot_list <- list()
    for (donor in good_donors) {
        plot_list[[donor]] <- spot_plot(
                spe, sample_id = donor, var_name = paste0("precast_k", k),
                is_discrete = TRUE
            ) +
                #   For some reason, spot_plot(colors = k_colors[[k]]) isn't
                #   working...
                scale_fill_manual(values = k_colors[[k]]) +
                #   Increase size of colored dots in legend
                guides(fill = guide_legend(override.aes = list(size = 5)))
    }

    pdf(file.path(plot_dir, sprintf("k%s_select_samples.pdf", k)))
    print(plot_list)
    dev.off()
}

session_info()
