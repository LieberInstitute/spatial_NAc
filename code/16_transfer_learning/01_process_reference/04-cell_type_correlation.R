####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)

spec <- matrix(
    c( "data", "d", 1, "character", "Specify input snRNA-seq dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)
opt <- list()
opt$data <- "human_NAc"
print(opt$data)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "preliminary_analysis", opt$data)
res_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning", "01_process_reference", "RCppML", opt$data)
dir.create(res_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

sce <- readRDS(file = file.path(dat_dir, "snRNA_seq_NAc.rds"))

x <- readRDS(file = file.path(res_dir,paste0("nmf_results.rds")))

####onehot encode cell type
if(opt$data == "human_NAc"){
    data<-as.data.frame(sce$CellType)
}else{
    data <- as.data.frame(sce$Combo_CellType)
}

colnames(data)<-'CellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ CellType, length)
rownames(onehot_cellType)<- onehot_cellType[,1]
onehot_cellType[ ,1] <- NULL
onehot_cellType <- onehot_cellType[match(rownames(t(x@h)) , rownames(onehot_cellType)), ]

###correlate with nmf patterns
pdf(file.path(plot_dir, paste0("nmf_cellType_correlation_heatmap.pdf")), width = 8, height = 10)
pheatmap(cor(t(x@h),onehot_cellType), fontsize_row = 9)
dev.off()

# aggregate NMF patterns

# create dataframe
if(opt$data == "human_NAc"){
    data <- data.frame(sce@meta.data, t(x@h))
    # aggregate NMF patterns across cell types
    # grep "NMF" to get all NMF patterns
    agg_data <- aggregate(data[,grep("nmf", colnames(data))],
                      by=list(data$CellType),
                      FUN=mean)
    # move Group.1 to row names, then drop
    rownames(agg_data) <- agg_data$Group.1
    agg_data <- agg_data[,-1]
}

if(opt$data == "rat_case_control_acute" | opt$data == "rat_case_control_repeated"){
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


pdf(file.path(plot_dir, paste0("nmf_cellType_correlation_aggregated_heatmap.pdf")), width = 10, height = 6)
p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column",
               fontsize_col = 9
)
dev.off()