####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(sessioninfo)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)
library(Seurat)

spec <- matrix(
    c("data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
#opt <- list()
#opt$data <- "human_NAc"
print(opt$data)

##load nmf patterns
x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', opt$data, paste0('nmf_results.rds')))

##load spe
raw_in_path <- here("processed-data", "05_harmony_BayesSpace", "02-compute_QC_metrics", "spe_with_QC_metrics_hdf5")
spe <- loadHDF5SummarizedExperiment(raw_in_path)
spe$low_umi <- spe$sum_umi < 250
spe$low_gene <- spe$sum_gene < 250
spe$low_gene_edge_spot <- spe$low_umi & spe$edge_distance < 6
spe <- spe[ ,spe$low_gene_edge_spot == "FALSE"]
spe <- spe[ ,spe$local_outliers == "FALSE"]
exprLogic <- counts(spe) > 0
spe$sample_id <- as.factor(spe$sample_id)
nSpots_by_donor <- lapply(levels(spe$sample_id), function(iSample){
        rowSums(exprLogic[ ,spe$sample_id == iSample])
})
nSpots_by_donor <- do.call(cbind, nSpots_by_donor)
colnames(nSpots_by_donor) <- levels(spe$sample_id)
select.genes <-  rownames(nSpots_by_donor)[rowSums(nSpots_by_donor == 0) == 0]
      # Only select those genes which have non-zero expression in atleast 1 spot in each slide
spe <- spe[rownames(spe) %in% select.genes, ]
counts <- as.matrix(counts(spe))
metadata <- colData(spe)
sp_obj <- CreateSeuratObject(counts = counts, project = "NAc_10x_visium")
metadata <- metadata[rownames(metadata) %in% colnames(sp_obj), ]
metadata <- metadata[match(colnames(sp_obj), rownames(metadata)), ]
sp_obj$donor <- metadata$sample_id
sp_obj$capture_area <- metadata$sample_id_original
sp_obj$slide_num <- metadata$slide_num
sp_obj$array_num <- metadata$array_num
sp_obj$age <- metadata$Age
sp_obj$sex <- metadata$Sex
sp_obj$exclude_overlapping <- metadata$exclude_overlapping
sp_obj$overlap_slide <- metadata$overlap_slide

sp_obj <- sp_obj %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)

logcounts(spe) <- sp_obj[["RNA"]]$data
#spe <- computeLibraryFactors(spe)
#spe <- logNormCounts(spe)
# Plot dir
plot_dir <- here("plots", "16_transfer_learning","02_target_projections", opt$data)
res_dir <- here("processed-data", "16_transfer_learning","02_target_projections", opt$data)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

##projection
set.seed(1029)
# extract patterns
patterns <- t(x@h)
colnames(patterns) <- paste("NMF", 1:dim(patterns)[2], sep = "_")

# extract loadings
loadings <- x@w

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated" | opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
    refDir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis")
    orthologs_df <- readRDS(file.path(refDir, opt$data, "orthologs_df.rds"))
    orthologs_df <- orthologs_df[match(rownames(loadings), orthologs_df$rat_genes), ]
    rownames(loadings) <- orthologs_df$human_genes
}
# ====== project loadings to spatial data =======
#rownames(spe) <-rowData(spe)$gene_name
if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"| opt$data == "rat_case_control_morphine_acute" | opt$data == "rat_case_control_morphine_repeated"){
    rownames(spe) <- make.names(rowData(spe)$gene_name, unique = TRUE)
}
common_genes <- intersect(rownames(loadings), rownames(spe))
loadings <- loadings[rownames(loadings) %in% common_genes,]
spe <- spe[rownames(spe) %in% common_genes,]
loadings <- loadings[match(rownames(spe),rownames(loadings)),]

logcounts <- logcounts(spe)

#loadings <- as(loadings, "dgCMatrix")

proj <- project(w=loadings, data=logcounts)

# remove rowSums == 0
proj1 <- proj[rowSums(proj) == 0, ,drop = FALSE]
proj2 <- proj[rowSums(proj) != 0, ,drop = FALSE]

proj2 <- apply(proj2,1,function(x){x/sum(x)})
proj1 <- t(proj1)

proj_final <- cbind(proj2, proj1)
proj_final <- proj_final[ ,match(rownames(proj), colnames(proj_final))]

colData(spe) <- cbind(colData(spe),proj_final)

saveRDS(spe, file = file.path(res_dir, paste0("spe_NMF.rds")))

################################################################

session_info()