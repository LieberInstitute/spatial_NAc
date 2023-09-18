library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)

source(here('code', '07_spot_deconvo', 'shared_functions.R'))

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds')
plot_dir = here('plots', '05_harmony_BayesSpace', 'test_plots')

donor = "Br2720" # "Br6471"
discrete_var = "10x_kmeans_7_clusters"
cont_var = "sum_gene"

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)

for (donor in unique(spe$sample_id)[unique(spe$sample_id) != "Br6432"]) {
    p = spot_plot(
        spe, donor, "", discrete_var, include_legend = FALSE, is_discrete = TRUE
    )
    pdf(file.path(plot_dir, sprintf('vis_merged_discrete_%s.pdf', donor)))
    print(p)
    dev.off()

    p = spot_plot(
        spe, donor, "", cont_var, include_legend = FALSE, is_discrete = FALSE
    )
    pdf(file.path(plot_dir, sprintf('vis_merged_continuous_%s.pdf', donor)))
    print(p)
    dev.off()
}
