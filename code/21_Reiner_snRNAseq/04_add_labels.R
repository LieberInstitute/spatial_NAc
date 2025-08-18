library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(scry)
library(here)

sce <- readRDS(file = here("processed-data","21_Reiner_snRNAseq","sce_clean.Rds"))

# Read in the meta data processed by the Reiner lab
annotated_df <- read.csv(here("processed-data", "21_Reiner_snRNAseq", "clean_metadata_w_labels.csv"))
rownames(annotated_df) <- annotated_df$X
annotated_df$X <- NULL

sce <- sce[ ,colnames(sce) %in% rownames(annotated_df)]
annotated_df <- annotated_df[match(colnames(sce), rownames(annotated_df)), ]
sce$group <- annotated_df$group
sce$subject <- annotated_df$subject
sce$treatment <- annotated_df$treatment
sce$cluster <- annotated_df$cluster

saveRDS(sce, file = here("processed-data","21_Reiner_snRNAseq","sce_clean_w_labels.Rds"))