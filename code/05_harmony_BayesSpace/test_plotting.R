library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
source(here('code', '05_harmony_BayesSpace', 'plotting_functions.R'))

spe_path = here('processed-data', '05_harmony_BayesSpace', 'spe.rds')
plot_dir = here('plots', '05_harmony_BayesSpace')

sample_id = "Br2720"
discrete_var = "X10x_kmeans_7_clusters"
cont_var = "sum_gene"

spe = readRDS(spe_path)

# colnames(spe) = spe$key
# spe$pxl_row_in_fullres_original <- spe$pxl_row_in_fullres
# spe$pxl_col_in_fullres_original <- spe$pxl_col_in_fullres
# spe$pxl_row_in_fullres <- spe$pxl_col_in_fullres <- NULL

p = vis_merged(spe, sampleid = sample_id, coldatavar = discrete_var)
pdf(file.path(plot_dir, 'vis_merged_discrete.pdf'))
print(p)
dev.off()

p = vis_merged(spe, sampleid = sample_id, coldatavar = cont_var)
pdf(file.path(plot_dir, 'vis_merged_continuous.pdf'))
print(p)
dev.off()
