######Load libraries
library(Seurat)
library(ggplot2)
library(Libra)
library(dplyr)
library(ggplot2)
library(stringr)
library(Matrix.utils)
library(ComplexUpset)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(here)

dat_dir <- here("processed-data", "12_snRNA")
res_dir <- here("processed-data", "17_cross_species_comparison", "rat_case_control")
plot_dir <- here("plots", "17_cross_species_comparison", "rat_case_control")
NAc_Combo <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))

#Create a sample id column that is Dataset_Sex_Stim
NAc_Combo$sample.id <- as.factor(paste(NAc_Combo$Dataset,NAc_Combo$Sex_Stim,sep = "_"))

#Get cell and sample metrics for aggregation
NAc_Combo@meta.data$CellType <- as.factor(NAc_Combo@meta.data$Combo_CellType)
cell_names <- purrr::set_names(levels(NAc_Combo@meta.data$CellType))

#Number of clusters
cluster_total <- length(cell_names)
# Named vector of sample names
NAc_Combo$Sex_Stim <- factor(NAc_Combo$Sex_Stim)
sample_ids <- purrr::set_names(levels(NAc_Combo@meta.data$Sex_Stim))
sample_ids
# Total number of samples 
sample_total <- length(sample_ids)
sample_total
#Figure out how many cells in each main group
table(NAc_Combo@meta.data$Sex_Stim)

##########Count aggregation to sample level
groups <- NAc_Combo@meta.data[, c("CellType", "sample.id")]

# Aggregate across cluster-sample groups (raw counts, thus counts slot)
count_aggr <- Matrix.utils::aggregate.Matrix(t(NAc_Combo@assays$RNA@counts), 
                                             groupings = groups, fun = "sum") 
dim(count_aggr)

#Transpose count_aggr
count_aggr_t <- t(count_aggr)

#change the count_aggr_t
colnames(count_aggr_t) <- gsub(x = colnames(count_aggr_t),pattern = "-",replacement = ".")

#Create a metadata dataframe
metadata <- data.frame(cluster_id = sub("_.*","", colnames(count_aggr_t)),
                       sample_id  = colnames(count_aggr_t),
                       Sex.Stim  = sub("^[^_]*_","", colnames(count_aggr_t)))
metadata$dataset <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",1)))
metadata$Sex     <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",2)))
metadata$Stim    <- as.factor(as.character(lapply(strsplit(metadata$Sex.Stim,"_"),"[",3)))
rownames(metadata) <- metadata$sample_id

###Calculate DEGs 
counts_mat <- NAc_Combo@assays$RNA@counts
Idents(NAc_Combo) <- gsub(x = Idents(NAc_Combo),pattern = "-",replacement = ".")
for(i in unique(metadata$cluster_id)){
  ############Make the DESEq object############
  metadata$CellType <- ifelse(metadata$cluster_id == i,
                              i,
                              "Other")
  dds <- DESeqDataSetFromMatrix(count_aggr_t, 
                                colData = metadata, 
                                design = ~ dataset + Stim + CellType)
  dds$CellType <- relevel(dds$CellType, ref = "Other")
  keep <- rowMeans(counts(dds)) > 5
  dds <- dds[keep,]
  print(paste("dds made for",i))
  ############Quality control############
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  #PCA for major celltype
  data <- plotPCA(vsd, intgroup = c("CellType"), returnData = TRUE)
  data$cluster_id <- as.character(lapply(strsplit(x = data$name,split = "_"),"[",1))
  percentVar <- round(100 * attr(data, "percentVar"))
  pca.plot.CellType <- ggplot(data, aes(PC1, PC2, color = CellType)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  #PCA for cluster_id
  pca.plot.cluster_id <- ggplot(data, aes(PC1, PC2, color = cluster_id)) +
    geom_point(size = 7) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme(text = element_text(size = 20)) +
    theme_bw(base_size = 16)
  ggsave(filename = paste0(plot_dir,"/", i, ".pdf"),
         plot     = cowplot::plot_grid(plotlist = list(pca.plot.CellType,pca.plot.cluster_id),ncol = 2),
         width = 15,
         height = 7)
  print(paste("PCA plots made for",i))
  ############Get results############
  dds <- DESeq(dds, test="LRT", reduced = ~ dataset + Stim)
  res <- as.data.frame(results(dds))
  res$GeneName <- rownames(res)
  cells.1    <- WhichCells(NAc_Combo,idents = i)
  cells.2    <- setdiff(x = row.names(NAc_Combo@meta.data),y=cells.1)
  res$Pct_Expressing <- NA
  res$Pct_CellType_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.1,drop = FALSE] > 0) / length(cells.1))*100
  res$Pct_Other_Expressing <- (rowSums(x = counts_mat[res$GeneName, cells.2,drop = FALSE] > 0) / length(cells.2))*100
  write.table(file      = paste0(res_dir, "/",i,".txt"),
              x         = res,
              sep       = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote     = FALSE)
  print(paste("DEG lists written for",i))
  rm(dds,keep,vsd,data,percentVar,pca.plot.CellType,pca.plot.cluster_id,res)
}

sessionInfo()