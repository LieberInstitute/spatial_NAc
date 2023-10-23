library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(Seurat)
library(SeuratData)

spe_in = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered.rds'
)
spot_diameter = 55e-6

#   Given a SpatialExperiment 'spe', convert to a Seurat object and return
spe_to_seurat = function(spe) {
    seur = as.Seurat(spe)

    return(seur)
}

spe = readRDS(spe_in)
seur_ex = LoadData("stxBrain", type = "anterior1") # InstallData("stxBrain")

this_sample = unique(spe$donor)[2]
spe_small = spe[1:100, spe$donor == this_sample]
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

class(seur_ex@images[[this_sample]]@image)
dim(seur_ex@images[[this_sample]]@image)
#   Values in [0, 1]

a = t(col2rgb(imgRaster(spe_small))) # makes vertical stripes
a = t(col2rgb(t(imgRaster(spe_small)))) # makes little repeated shapes?
a = col2rgb(imgRaster(spe_small)) # vertical stripes with pink on left and blue on right
a = col2rgb(t(imgRaster(spe_small))) # RGB diagonal stripes

this_img = array(
    a,
    dim = c(dim(imgRaster(spe_small)), 3)
) / 256

seur_spe@images[[this_sample]] = Seurat:::VisiumV1(
    image = this_img,
    scale.factors = scalefactors(
        spot = NA, fiducial = NA, hires = NA, lowres = scaleFactors(spe_small) / 3
    ),
    coordinates = coords,
    spot.radius = spot_diameter / scaleFactors(spe_small),
    assay = "Spatial",
    key = paste0(this_sample, '_')
)
names(seur_spe@images) = this_sample # why is this necessary?

# SpatialFeaturePlot(seur_ex, features = c("Hpca", "Ttr"))
SpatialFeaturePlot(seur_spe, features = rownames(seur_spe)[1:2])#, images = this_sample)
