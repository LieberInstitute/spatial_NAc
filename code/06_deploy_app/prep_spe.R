library("SpatialExperiment")
library("lobstr")

spe_path_in = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
spe_path_out = here(
    'processed-data', '06_deploy_app', 'spe_shiny.rds'
)

dir.create(dirname(spe_out_path), showWarnings = FALSE)

spe = readRDS(spe_path)

#   Drop counts assay and check size in memory
assays(spe)$counts = NULL
size_string = paste0(round(obj_size(spe) / 1e9, 1), 'GB')
paste('Object size:', size_string)

saveRDS(spe, spe_path_out)
