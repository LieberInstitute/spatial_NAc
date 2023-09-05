library("SpatialExperiment")
library("lobstr")

spe_path_in = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
spe_path_out <- here(
    "code", "06_deploy_app", "spe_shiny.rds"
)

dir.create(dirname(spe_path_out), showWarnings = FALSE, recursive = TRUE)

spe <- readRDS(spe_path_in)
spe
# class: SpatialExperiment
# dim: 30854 135597
# metadata(0):
# assays(2): counts logcounts
# rownames(30854): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
#   ENSG00000277196
# rowData names(7): source type ... gene_type gene_search
# colnames(135597): AAACAACGAATAGTTC-1_V11D01-384_B1
#   AAACAAGTATCTCCCA-1_V11D01-384_B1 ... TTGTTTGTATTACACG-1_V13M06-378_D1
#   TTGTTTGTGTAAATTC-1_V13M06-378_D1
# colData names(41): sample_id in_tissue ... scran_quick_cluster
#   sizeFactor
# reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# imgData names(4): sample_id image_id data scaleFactor
obj_size(spe)
# 6.99 GB

#   Drop counts assay and check size in memory
assays(spe)$counts = NULL
size_string = paste0(round(obj_size(spe) / 1e9, 1), 'GB')
paste('Object size:', size_string)

saveRDS(spe, spe_path_out)
