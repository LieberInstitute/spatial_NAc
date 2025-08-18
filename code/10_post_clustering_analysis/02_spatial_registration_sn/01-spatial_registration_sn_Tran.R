library(SpatialExperiment)
library(spatialLIBD)
library(jaffelab)
library(here)
library(sessioninfo)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(fastTopics)
library(HDF5Array)
library(tidyverse)
library(getopt)
library(edgeR)
library(scater)
library(scran)
library(DeconvoBuddies)
library(dplyr)
library(gridExtra)
library(ggforce)
library(pheatmap)

# Load the data
load(here("processed-data", "SCE_NAc-n8_tran-etal.rda"))
sce <- sce.nac.tran
rm(sce.nac.tran)

sce <- sce[ ,!grepl("drop", sce$cellType)]
sce$cellType <- as.character(sce$cellType)

subset_neurons <- FALSE
if(subset_neurons){
        select_cell_types <- c("Inhib_A", "Inhib_B", "Inhib_C", "Inhib_D", "Inhib_E", "MSN.D1_A", 
        "MSN.D1_B", "MSN.D1_C", "MSN.D1_D", "MSN.D1_E", "MSN.D1_F", "MSN.D2_A", "MSN.D2_B", "MSN.D2_C", 
        "MSN.D2_D")
        sce <- sce[ ,sce$cellType %in% select_cell_types]
}
sce$cellType <- factor(sce$cellType)

res_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
plot_dir <- here("plots", "10_post_clustering_analysis", "02_spatial_registration_sn")

counts <- counts(sce)
cellType_abundance <- data.frame(table(sce$cellType))
colnames(cellType_abundance) <- c("Cell_type", "Ncells")
if(subset_neurons){
        pdf(file.path(plot_dir, "cellType_abundance_Tran_neurons.pdf"), width = 12, height = 5)
        ggplot(cellType_abundance, aes(x = Cell_type, y = Ncells, fill = Cell_type)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cell Type") + ylab("Number of cells") + coord_flip() + theme(legend.position = "none") 
        dev.off()
}else{
        pdf(file.path(plot_dir, "cellType_abundance_Tran.pdf"), width = 12, height = 5)
        ggplot(cellType_abundance, aes(x = Cell_type, y = Ncells, fill = Cell_type)) + geom_bar(stat = "identity") + theme_classic() + xlab("Cell Type") + ylab("Number of cells") + coord_flip() + theme(legend.position = "none") 
        dev.off()
}

## Factor categorical variables used as covariates
colData(sce)$sex <- factor(colData(sce)$sex, levels = c("F", "M"))
## Use all unique ensembl IDs as rownames
rownames(sce) <- rowData(sce)$gene_id

markers_1vALL_enrich <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "cellType", 
                                          mod = "~donor")
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")
if(subset_neurons){
        save(markers_1vALL_df,file = file.path(res_dir, "markers_1vAll_CellType_Final_tran_neurons.rda"))
}else{
        save(markers_1vALL_df,file = file.path(res_dir, "markers_1vAll_CellType_Final_tran.rda"))
}


# Check the drivers of variation in pseudobulk
sce_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        cluster = colData(sce)$cellType,
        sample_id = colData(sce)$donor
    ))
sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]
colData(sce_pseudo)<-colData(sce_pseudo)[,c('region', 'donor','sex','processBatch','protocol','sequencer',
                                            'cellType','ncells', 'cluster', 'sample_id')]

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_sample_id <- filterByExpr(sce_pseudo, group = sce_pseudo$sample_id)
rowData(sce_pseudo)$high_expr_group_cellType <- filterByExpr(sce_pseudo, group = sce_pseudo$cellType)

sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_sample_id, ]
x <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)
stopifnot(identical(rownames(x), rownames(sce_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(sce_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
logcounts(sce_pseudo) <- x
sce_pseudo <- scuttle::addPerCellQC(
    sce_pseudo,
    subsets = list(Mito = which(seqnames(sce_pseudo) == "chrMT"))
)
sce_pseudo$det_out<-as.logical(isOutlier(sce_pseudo$detected,type='lower',nmads=3))
sce_pseudo<-sce_pseudo[,sce_pseudo$det_out==F]

pca <- prcomp(t(assays(sce_pseudo)$logcounts))

message(Sys.time(), " % of variance explained for the top 20 PCs:")
metadata(sce_pseudo)
metadata(sce_pseudo) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
metadata(sce_pseudo)

pca_pseudo <- pca$x[, seq_len(50)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

jaffelab::getPcaVars(pca)[seq_len(50)]

if(subset_neurons){
        pdf(file = file.path(plot_dir ,"pseudobulk_PCA_sce_final_tran_neurons.pdf"), width = 14, height = 14)
plotPCA(sce_pseudo, colour_by = "cellType", ncomponents = 5, point_size = 1, point_alpha=1,label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "donor", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "sex", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "protocol", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "processBatch", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
dev.off()
}else{
        pdf(file = file.path(plot_dir ,"pseudobulk_PCA_sce_final_tran.pdf"), width = 14, height = 14)
plotPCA(sce_pseudo, colour_by = "cellType", ncomponents = 5, point_size = 1, point_alpha=1,label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "donor", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "sex", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "protocol", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "processBatch", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "subsets_Mito_percent", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "sum", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
plotPCA(sce_pseudo, colour_by = "detected", ncomponents = 5, point_size = 1, label_format = c("%s %02i", " (%i%%)"),
        percentVar = metadata(sce_pseudo)$PCA_var_explained)
dev.off()
}


## Run single cell registration to identify markers
sn_registration <- registration_wrapper(
    sce = sce,
    var_registration = "cellType",
    var_sample_id = "donor",
    covars = c("sex", "processBatch"),
    gene_ensembl = "gene_id",
    gene_name = "gene_name"
)

if(subset_neurons){
        saveRDS(sn_registration, file = file.path(res_dir, "sn_cellType_registration_tran_neurons.rds"))
}else{
        saveRDS(sn_registration, file = file.path(res_dir, "sn_cellType_registration_tran.rds"))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()