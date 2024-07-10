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
        "marker_genes", "m", 1, "logical", "Use only marker genes from single cell labels?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

# Read in spatial data
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5")
spe <- loadHDF5SummarizedExperiment(spe_dir)
sample_id <- levels(spe$sample_id)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe), unique = TRUE)
spe <- spe[ ,spe$sample_id == sample_id] 

# Specify location to save files and plots
processed_out <- here("processed-data","20_RCTD")
plots_out <- here("plots","20_RCTD")

# Create spatial RNA object, which requires, counts, coords, and nUMI
counts = assays(spe)$counts
coords = as.data.frame(spatialCoords(spe))
rownames(coords) <- colnames(counts)
colnames(coords) = c("xcoord","ycoord")
nUMI <- colSums(counts) 
puck <- SpatialRNA(coords, counts, nUMI)

# Plot diagnostics and save
barcodes <- colnames(puck@counts)
dir.create(plots_out, showWarnings = FALSE)
dir.create(here(plots_out, sample_id))

pdf(here(plots_out, sample_id, "UMIcount.pdf"), width = 6, height = 6)
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
dev.off()

pdf(here(plots_out, sample_id, "Distribution_log_UMI.pdf"), width = 6, height = 6)
hist(log(puck@nUMI,2), xlab = "log2(nUMI)", main = "Distribution of log2(nUMI)")
dev.off()

if(opt$marker_genes){
    reference = readRDS(here(processed_out,'SCRef_markers.rds'))
}else{
    reference = readRDS(here(processed_out,'SCRef.rds'))
}

myRCTD <- create.RCTD(puck, reference, max_cores = 1, MAX_MULTI_TYPES = 5, UMI_min = 2)

dir.create(here(processed_out, sample_id))

if(opt$marker_genes){
    saveRDS(myRCTD, here(processed_out,sample_id,'myRCTD_markers.rds'))
}else{
    saveRDS(myRCTD, here(processed_out,sample_id,'myRCTD.rds'))
}

