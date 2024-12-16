
####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

spec <- matrix(
    c(
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs", 
        "data", "d", 1, "character", "Specify input snRNA-seq dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
opt <- list()
opt$gene_selection_strategy <- "all_genes"
opt$data <- "rat_case_control"
print(opt$gene_selection_strategy)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
  }else{
        stop("Invalid input data set")
  }
}

x <- readRDS(file = file.path(res_dir,paste0("nmf_results_",opt$gene_selection_strategy, ".rds")))

####onehot encode cell type
if(opt$data == "human_NAc"){
    data<-as.data.frame(sce$CellType.Final)
}else{
    data <- as.data.frame(sce$Combo_CellType)
}

colnames(data)<-'CellType_final'
onehot_cellType_final <-  dcast(data = data, rownames(data) ~ CellType_final, length)
rownames(onehot_cellType_final)<- onehot_cellType_final[,1]
onehot_cellType_final[ ,1] <- NULL
onehot_cellType_final <- onehot_cellType_final[match(rownames(t(x@h)) , rownames(onehot_cellType_final)), ]

###correlate with nmf patterns
pdf(file.path(plot_dir, paste0("nmf_cellType_correlation_heatmap_", opt$gene_selection_strategy,".pdf")))
pheatmap(cor(t(x@h),onehot_cellType_final), fontsize_row = 9)
dev.off()

# aggregate NMF patterns

# create dataframe
if(opt$data == "human_NAc"){
    data <- data.frame(colData(sce), t(x@h))
    # aggregate NMF patterns across cell types
    # grep "NMF" to get all NMF patterns
    agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$CellType.Final),
                      FUN=mean)
    # move Group.1 to row names, then drop
    rownames(agg_data) <- agg_data$Group.1
    agg_data <- agg_data[,-1]
}

if(opt$data == "rat_case_control"){
    data <- data.frame(sce@meta.data, t(x@h))
    # aggregate NMF patterns across cell types
    # grep "NMF" to get all NMF patterns
    agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$Combo_CellType),
                      FUN=mean)
    # move Group.1 to row names, then drop
    rownames(agg_data) <- agg_data$Group.1
    agg_data <- agg_data[,-1]
}


pdf(file.path(plot_dir, paste0("nmf_cellType_correlation_aggregated_heatmap_", opt$gene_selection_strategy,".pdf")))
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 9
)
dev.off()