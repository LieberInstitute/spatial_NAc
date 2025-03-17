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

spec <- matrix(
    c(  "data", "d", 1, "character", "Specify the dataset to be used?",
        "gene_selection_strategy", "g", 1, "character", "Choose all genes, or highly deviant genes based on snRNA-seq, or nnSVGs"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

opt <- list()
opt$data <- "human_NAc"

print(opt$data)
res_dir <- here::here("processed-data", "16_transfer_learning","01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning","01_process_reference", "RCppML", opt$data)

go <- readRDS(file=file.path(res_dir,'go_analysis.rds'))
names(go) <- paste0("NMF_", c(1:length(go)))

# Subset only to NMFs that have some pathways that are significantly over-represented
length_enr_res <- lapply(go, function(igo){
    length(igo$Description)
})
length_enr_res <- unlist(length_enr_res)
go <- go[length_enr_res > 0]

plot_list <- list()
for(i in c(1:length(go))){
    igo <- go[[i]]
    plot_list[[i]] <- dotplot(igo, showCategory=30) + ggtitle(gsub("_", " ", names(go)[i]))
}
pdf(file.path(plot_dir, "nmf_gene_ora_analysis.pdf"), height = 18, width = 9)
print(plot_list)
dev.off()
