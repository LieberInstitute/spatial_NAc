library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(Seurat)
library(SeuratData)

spe_in = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)

#   Given a SpatialExperiment 'spe', convert to a Seurat object and return
spe_to_seurat = function(spe) {
    seur = as.Seurat(spe)

    return(seur)
}

spe = readRDS(spe_in)
seur_ex = LoadData("stxBrain", type = "anterior1") # InstallData("stxBrain")

spe_small = spe[1:100, spe$donor == unique(spe$donor)[2]]
seur_spe = as.Seurat(spe_small)

#   Convert 'spatialCoords' slot into a Seurat-compatible equivalent
coords = colData(spe_small)[
    ,
    c(
        'in_tissue', 'array_row_transformed', 'array_col_transformed',
        'pxl_row_in_fullres_transformed', 'pxl_col_in_fullres_transformed'
    )
]
colnames(coords) = c('tissue', 'row', 'col', 'imagerow', 'imagecol')
coords$tissue = as.integer(coords$tissue)
coords = as.data.frame(coords)

Seurat:::VisiumV1(
    image = as.matrix(imgRaster(spe_small)),
    scale.factors = scalefactors(
        spot = NA, fiducial = NA, hires = NA, scaleFactors(spe_small)
    ),
    coordinates = coords,
    spot.radius = 27.5e-6
)

SpatialFeaturePlot(seur_ex, features = c("Hpca", "Ttr"))
