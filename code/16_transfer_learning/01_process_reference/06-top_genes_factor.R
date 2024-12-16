library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
library(here)
library(scuttle)
library(dplyr)
library(getopt)
library(Matrix)
library(HDF5Array)
set.seed(1234)

spec <- matrix(
    c(
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs", 
        "data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

#opt <- list()
#opt$gene_selection_strategy <- "all_genes"
#opt$data <- "human_NAc"
print(opt$gene_selection_strategy)
print(opt$data)

# Read data based on whether we are processing the human or rat NAc snRNA seq data
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

x <- readRDS(file = file.path(res_dir,paste0("nmf_results_", opt$gene_selection_strategy, ".rds")))

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

geneData <- rowData(sce)
geneData <- geneData[geneData$gene_id %in% rownames(sce), ]
geneData <- geneData[match(rownames(x@w), geneData$gene_id), ]
rownames(x@w) <- geneData$gene_name
# get top 50 genes
top50 <- top_genes(x@w, 50)
write.csv(top50, file = file.path(res_dir, "top50_genes.csv"))

rownames(sce) <- rowData(sce)$gene_name
sce <- sce[rownames(sce) %in% rownames(x@w), ]
expr <- logcounts(sce)


corr_mat <- matrix(NA, nrow = dim(sce)[1], ncol = dim(x@h)[1])
for(i in c(1:dim(x@h)[1])){
  cat(i, "\n")
  corr_mat[ ,i] <- apply(expr, 1, function(igene) cor(igene, x@h[i, ]))
}
rownames(corr_mat) <- rownames(expr)
colnames(corr_mat) <- rownames(x@h)

saveRDS(corr_mat, file.path(res_dir, "gene_corr.rds"))
corr_mat <- reshape2::melt(corr_mat)
corr_mat <- corr_mat[abs(corr_mat$value) > 0.2, ]
colnames(corr_mat) <- c("gene", "factor", "correlation")

write.csv(corr_mat, file.path(res_dir, paste0("gene_corr_thresholded_", opt$gene_selection_strategy, ".csv")), row.names= FALSE, quote = FALSE)