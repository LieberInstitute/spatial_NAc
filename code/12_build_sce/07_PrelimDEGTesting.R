#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(scran)
library(here)

##############
#Load sce
print("Loading SCE")
Sys.time()
sce <- readRDS(here("processed-data","12_snRNA","sce_clustered_log.Rds"))

dim(sce)

stopifnot(identical(rownames(colData(sce)),colnames(sce)))

sce
###############
print("Starting DEG testing")
mod <- with(colData(sce), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(sce, 
                                groups=sce$k_10_louvain_1,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})

#Add gene info to each list.
for(i in names(markers_pairwise)){
    markers_pairwise[[i]] <- as.data.frame(markers_pairwise[[i]])
    markers_pairwise[[i]]$gene_id <- row.names(markers_pairwise[[i]])
    markers_pairwise[[i]] <- dplyr::left_join(x  =  markers_pairwise[[i]],
                                              y  =  as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                              by = "gene_id")
}


save(markers_pairwise,file = here("processed-data","12_snRNA","markers_pairwise_list_k_10_louvain_1.rda"))

########
#1vALL DEG testing
markers_1vALL_enrich <- findMarkers_1vAll(sce, 
                                          assay_name = "logcounts", 
                                          cellType_col = "k_10_louvain_1", 
                                          mod = "~Sample")

#Add symbol information to the table
#First change the ensembl gene id column to have same name as what is in rowData(sce)
colnames(markers_1vALL_enrich)[1] <- "gene_id"
markers_1vALL_df <- dplyr::left_join(x = as.data.frame(markers_1vALL_enrich),
                                     y = as.data.frame(rowData(sce)[,c("gene_id","gene_name")]),
                                     by = "gene_id")

#save the dataframe. 
save(markers_1vALL_df,file = here("processed-data","markers_1vAll_k_10_louvain_1.rda"))


print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
