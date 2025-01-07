library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

#load("sce.Rdata", verbose = TRUE)
sce <- readRDS(file = "sce_NAc_app.Rds")
#load("sce_for_iSEE_LS.rda", verbose = TRUE)

#Change the rownames frome ensembl id to gene_name
rownames(sce) <-  uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)

#Source
source("initial.R", print.eval = TRUE)

#Increase minimum number of colors to the # of clusters within CellType.Final
sce <- registerAppOptions(sce, color.maxlevels = length(unique(sce$CellType.Final)))

#Deploy app
iSEE(
    sce,
    appTitle = "NAc snRNAseq data",
    initial = initial
)
