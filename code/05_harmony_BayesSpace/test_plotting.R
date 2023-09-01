library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
source(here('code', '05_harmony_BayesSpace', 'plotting_functions.R'))

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds')
plot_dir = here('plots', '05_harmony_BayesSpace', 'test_plots')

sample_id = "Br2720"
discrete_var = "10x_kmeans_7_clusters"
cont_var = "sum_gene"

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)

for (donor in unique(spe$sample_id)[unique(spe$sample_id) != "Br6432"]) {
    p = vis_merged(spe, sampleid = donor, coldatavar = discrete_var)
    pdf(file.path(plot_dir, sprintf('vis_merged_discrete_%s.pdf', donor)))
    print(p)
    dev.off()

    p = vis_merged(spe, sampleid = donor, coldatavar = cont_var)
    pdf(file.path(plot_dir, sprintf('vis_merged_continuous_%s.pdf', donor)))
    print(p)
    dev.off()
}
