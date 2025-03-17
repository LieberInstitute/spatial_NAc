####code for correlating categorical technical variables with NMF patterns
library(here)
library(SpatialExperiment)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(dplyr)
library(RcppML)
library(getopt)
library(org.Hs.eg.db)
library(clusterProfiler)


source(here("code", "16_transfer_learning","01_process_reference", "utilities.R"))

spec <- matrix(
    c(  "data", "d", 1, "character", "Specify the dataset to be used?"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

opt <- list()
opt$data <- "human_NAc"
print(opt$data)

# Read data and create Seurat object
dat_dir <- here::here("processed-data", "12_snRNA")
res_dir <- here::here("processed-data", "16_transfer_learning","01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning","01_process_reference", "RCppML", opt$data)

if(opt$data == "human_NAc"){
  sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))
}else{
  if(opt$data == "rat_case_control"){
    sce <- readRDS(file = file.path(dat_dir, "NAc_Combo_Integrated.RDS"))
  }else{
        stop("Invalid input data set")
  }
}

x <- readRDS(file = file.path(res_dir,paste0("nmf_results.rds")))

## Set up marker gene detection
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# Remove mitochondrial genes
mito <- rownames(sce)[which(seqnames(sce) == "chrM")]
loads <- loads[!rownames(loads) %in% mito,]
## now get markers
markers <- patternMarkers(loads,x@h[rownames(x@h) %in% colnames(loads),],'all',1,100)

genes <- markers$PatternMarkers
names(genes) <- colnames(loads)

go<-list()
for(i in 1:length(genes)){
cat(i, "\n")
go[[i]] <- enrichGO(gene  = genes[[i]],
                    universe      = rownames(loads),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1,
                    readable      = TRUE,
                    keyType = 'ENSEMBL')
}

### save go
saveRDS(go,file=file.path(res_dir,'go_analysis.rds'))
