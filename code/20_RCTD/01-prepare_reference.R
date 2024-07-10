library(spacexr)
library(Matrix)
library(SingleCellExperiment)
library(here)
library(scran)
library(scater)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(HDF5Array)
library(ggplot2)
library(patchwork)
library(getopt)

spec <- matrix(
    c(
        "marker_genes", "m", 1, "logical", "Use only marker genes from single cell labels?",
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

# Read and process snRNA-seq data
dat_dir <- here::here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
colData(sce)$CellType.Final[colData(sce)$CellType.Final == "T-Cell"] <- "T_cell"

# Find marker genes if we need to subset to only include genes that distinguish between cell-types
# Rename the genes using symbols instead of ensembl IDs
rownames(sce) <- rowData(sce)$gene_name
rownames(sce) <- make.names(rownames(sce), unique = TRUE)
sce <- logNormCounts(sce)
# Get vector indicating which genes are neither ribosomal or mitochondrial
genes <- !grepl(pattern = "^RP[SL]|MT\\.", x = rownames(sce))

# Obtain marker genes for each cell identity
colLabels(sce) <- colData(sce)$CellType.Final
# Compute marker genes
mgs <- scoreMarkers(sce, subset.row = genes)

# Only keep genes that are relevant for each cell identity
mgs_fil <- lapply(names(mgs), function(i) {
    cat(i, "\n")
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.75
    x <- x[x$mean.AUC > 0.75, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
})
mgs_df <- do.call(rbind, mgs_fil)
marker_genes <- unique(mgs_df$gene)

if(opt$marker_genes){
    sce <- sce[rownames(sce) %in% marker_genes, ]
}else{
    sce <- sce[genes, ]
}
print(dim(sce))

# Extract data to build the reference
counts <- counts(sce)
meta_data <- colData(sce)
cell_types <- meta_data$CellType.Final
names(cell_types) <- colnames(sce)
cell_types <- factor(cell_types, levels = unique(cell_types))
nUMI <- meta_data$sum
names(nUMI) <- colnames(sce)

# Create the reference object
refdir <- here::here("processed-data", "20_RCTD")
reference <- Reference(counts, cell_types, nUMI, n_max_cells = dim(sce)[2])

print(dim(reference@counts))
table(reference@cell_types)

if(opt$marker_genes){
    saveRDS(reference, file.path(refdir,'SCRef_markers.rds'))
}else{
    saveRDS(reference, file.path(refdir,'SCRef.rds'))
}
