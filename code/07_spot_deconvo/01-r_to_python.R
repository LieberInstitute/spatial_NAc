#   Convert SpatialExperiment and snRNA-seq SingleCellExperiment to AnnDatas.
#   Save a copy of the modified SCE as an R object as well. This is a processing
#   step prior to running cell2location. Finally, run getMeanRatio2 to rank
#   genes as markers, and save the resulting object.

library("basilisk")
library("SingleCellExperiment")
library("SpatialExperiment")
library("DeconvoBuddies")
library("zellkonverter")
library("sessioninfo")
library("here")
library("tidyverse")

cell_type_var = "TODO"
donor_formula = "TODO"

#  Paths
sce_in <- "TODO"
spe_in <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered.rds"
)

sce_out <- here("processed-data", "07_spot_deconvo", "sce.h5ad")
sce_r_out <- here("processed-data", "07_spot_deconvo", "sce.rds")
spe_out <- here("processed-data", "07_spot_deconvo", "spe.h5ad")

# sample_IF_out <- here(
#     "processed-data", "spot_deconvo", "05-shared_utilities", "IF", "sample_ids.txt"
# )
# sample_nonIF_out <- here(
#     "processed-data", "spot_deconvo", "05-shared_utilities", "nonIF", "sample_ids.txt"
# )
marker_object_out <- here("processed-data", "07_spot_deconvo", "marker_stats.rds")

#  Make sure output directory exists
dir.create(dirname(spe_out), recursive = TRUE, showWarnings = FALSE)

###############################################################################
#  Functions
###############################################################################

write_anndata <- function(sce, out_path) {
    invisible(
        basiliskRun(
            fun = function(sce, filename) {
                library("zellkonverter")
                library("reticulate")

                # Convert SCE to AnnData:
                adata <- SCE2AnnData(sce)

                #  Write AnnData object to disk
                adata$write(filename = filename)

                return()
            },
            env = zellkonverterAnnDataEnv(),
            sce = sce,
            filename = out_path
        )
    )
}

###############################################################################
#   Main
###############################################################################

#   Load objects
load(sce_in, verbose = TRUE)
spe <- readRDS(spe_in)
gc()

#-------------------------------------------------------------------------------
#   Convert snRNA-seq and spatial R objects to AnnData python objects
#-------------------------------------------------------------------------------

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe)$spatial <- spatialCoords(spe)

#   Use Ensembl gene IDs for rownames (not gene symbol)
rownames(sce) <- rowData(sce)$gene_id

#   Save a copy of the filtered + slightly modified sce as an R object, and
#   convert all objects to Anndatas
saveRDS(sce, sce_r_out)

print("Converting objects to AnnDatas...")
write_anndata(sce, sce_out)
write_anndata(spe, spe_out)
gc()

#-------------------------------------------------------------------------------
#   Rank marker genes
#-------------------------------------------------------------------------------

#   We won't consider genes that aren't in both objects
keep <- rowData(sce)$gene_id %in% rowData(spe)$gene_id
sce <- sce[keep, ]

perc_keep <- 100 * (1 - length(which(keep)) / length(keep))
print(
    paste0(
        "Dropped ", round(perc_keep, 1), "% of potential marker genes ",
        "that were not present in the spatial data"
    )
)

print("Running getMeanRatio2 and findMarkers_1vAll to rank genes as markers...")
marker_stats <- get_mean_ratio2(
    sce,
    cellType_col = cell_type_var, assay_name = "logcounts"
)
marker_stats_1vall <- findMarkers_1vAll(
    sce,
    cellType_col = cell_type_var, assay_name = "logcounts",
    mod = donor_formula
)
marker_stats <- left_join(
    marker_stats, marker_stats_1vall,
    by = c("gene", "cellType.target")
)

saveRDS(marker_stats, marker_object_out)

session_info()
