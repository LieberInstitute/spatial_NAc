library(here)
library(tidyverse)
library(ggplot2)
library(spatialLIBD)
library(sessioninfo)
library(HDF5Array)

sample_info_path = here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv')
plot_dir = here('plots', 'VistoSeg')
spe_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)

dir.create(plot_dir, showWarnings = FALSE)

sample_info = read.csv(sample_info_path) |> as_tibble()

################################################################################
#   Read in outputs from VistoSeg countNuclei for all samples and concatenate
#   into one tibble
################################################################################

spe_nac = loadHDF5SummarizedExperiment(spe_dir)

counts_list = list()
for (i in 1:nrow(sample_info)) {
    counts_list[[i]] = read.csv(
                file.path(sample_info$spaceranger_dir[i], 'tissue_spot_counts.csv')
            ) |>
            as_tibble() |>
            mutate(sample_id = sample_info$sample_id[i])
}

#   We aren't interested in spots not in tissue (e.g. spots in holes in the
#   tissue), and don't want them throwing off median, etc calculations
counts = do.call(rbind, counts_list) |>
    mutate(key = paste(barcode, sample_id, sep = '_')) |>
    filter(key %in% spe_nac$key[spe_nac$in_tissue])

################################################################################
#   Exploratory plots checking the distribution of nuclei counts across spots
#   and samples
################################################################################

#   Individual boxplots, one per sample
p = ggplot(counts) +
    geom_boxplot(
        aes(x = sample_id, y = Nmask_dark_blue),
        outlier.shape = NA
    ) +
    coord_cartesian(
        ylim = c(0, 2 * max(boxplot.stats(counts$Nmask_dark_blue)$stats))
    )
pdf(file.path(plot_dir, 'nuclei_counts_by_sample.pdf'))
print(p)
dev.off()

#   Calculate an upper cutoff for outliers, including all spots in all samples
upper_cutoff = median(counts$Nmask_dark_blue) + 3 * sd(counts$Nmask_dark_blue)

p = counts |>
    #   Determine the proportion, by sample, of spots whose counts exceed the
    #   cutoff
    group_by(sample_id) |>
    summarize(
        prop_outliers = mean(Nmask_dark_blue > upper_cutoff)
    ) |>
    #   Label a sample 'Not in top 5' if it doesn't have one of the top-5
    #   outlier proportions
    mutate(
        sample_id_plot = ifelse(
            prop_outliers %in% sort(prop_outliers, decreasing = TRUE)[1:5],
            sample_id,
            'Not in top 5'
        )
    ) |>
    #   Violin plot of outlier proportion by sample
    ggplot(aes(x = 1, y = prop_outliers)) +
        geom_violin() +
        geom_jitter(aes(color = sample_id_plot)) +
        labs(color = 'Sample ID', x = '', y = 'Proportion of Outlier Spots')
pdf(file.path(plot_dir, 'prop_outlier_spots_by_sample.pdf'))
print(p)
dev.off()

#   Verify that distribution of nuclei counts is somewhat comparable to those
#   observed in DLPFC
spe_dlpfc = fetch_data(type = "spatialDLPFC_Visium")
stopifnot(all(spe_dlpfc$in_tissue))
message(
    sprintf(
        "Median number of nuclei per spot (NAc, DLPFC): %s, %s",
        median(counts$Nmask_dark_blue),
        median(spe_dlpfc$VistoSeg_count)
    )
)

message(
    sprintf(
        "3 SDs above median count for nuclei per spot (NAc, DLPFC): %s, %s",
        round(upper_cutoff, 1),
        round(median(spe_dlpfc$VistoSeg_count) + 3 * sd(spe_dlpfc$VistoSeg_count))
    )
)

session_info()
