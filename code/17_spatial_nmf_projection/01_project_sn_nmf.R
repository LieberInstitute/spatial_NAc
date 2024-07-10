####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
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

##load nmf patterns
NMF_results <- readRDS(file=here::here('processed-data','19_snRNAseq_NMF','RcppML','top2000_HVGs', 'nmf_model.rds'))

##load spe
spe_in_path <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)
spe <- loadHDF5SummarizedExperiment(spe_in_path)

# Plot dir
plot_dir <- here("plots", "17_spatial_nmf_projection")
res_dir <- here("processed-data", "17_spatial_nmf_projection")
# Add PRECAST results to the spe
precast_path <- here("processed-data", "10_precast", "nnSVG_precast", "PRECAST_k%s.csv")
result_list <- list()
for (K in 3:28) {
    result_list[[K]] <- sprintf(precast_path, K) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster)) |>
        mutate(k = K)
}
precast_results <- do.call(rbind, result_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "precast_k"
    )

precast_results <- data.frame(precast_results)
for(i in c(2:dim(precast_results)[2])){
    precast_results[ ,i] <- factor(as.character(precast_results[ ,i]), levels = c(1:(i+1)))
}

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

##projection
set.seed(1029)
#loadings <- NMF_results@w %*% Diagonal(x = NMF_results@d)
loadings <- NMF_results@w
#loadings <- as.matrix(loadings)
rownames(spe) <- rowData(spe)$gene_name
common.genes <- intersect(rownames(spe),rownames(loadings))
loadings <- loadings[rownames(loadings) %in% common.genes,]
spe2 <- spe[rownames(spe) %in% common.genes,]
loadings <- loadings[match(rownames(spe2),rownames(loadings)),]

# First learn a simple single step projection
proj <- project(w = loadings, data = logcounts(spe2), L1=0.00005)
proj <- t(proj)
# Adding simple projection to spe object
colData(spe) <- cbind(colData(spe), proj)
# For each nmf plot the number of spots with non-zero weight in each sample
pdf(file.path(plot_dir, "number_spots_with_nonzero_coefficients.pdf"), width = 9, height = 5)
p_list <- list()
for(i in c(1:dim(proj)[2])){
    spot_freq <- data.frame(table(spe$sample_id[colData(spe)[ ,paste0("nmf", i)] > 0]))
    p_list[[i]] <- ggplot(spot_freq, aes(x = Var1, y = Freq, fill = Var1))+ geom_bar(stat = "identity")+ theme_classic() + 
    xlab("Sample ID") + ylab("Number of spots with non-zero loadings") + guides(fill=guide_legend(title="Sample ID")) + ggtitle(paste0("NMF ", i))
    print(p_list[[i]])
}
dev.off()

spot_plot(spe, "Br8325", var_name = "ENSG00000197971", is_discrete = FALSE, spatial = TRUE, min = 0)
spot_plot(spe, "Br6522", var_name = "ENSG00000197971", is_discrete = FALSE, spatial = TRUE, min = 0)
spot_plot(spe, "Br3942", var_name = "ENSG00000197971", is_discrete = FALSE, spatial = TRUE, min = 0)

spot_plot(spe, "Br2720", var_name = "nmf5", is_discrete = FALSE, spatial = TRUE, min = 0)
spot_plot(spe, "Br6522", var_name = "nmf45", is_discrete = FALSE, spatial = TRUE, min = 0)
spot_plot(spe, "Br3942", var_name = "nmf45", is_discrete = FALSE, spatial = TRUE, min = 0)


df <- cbind(colData(spe), proj)
df <- data.frame(df)
ggplot(df, aes(x = nmf6, y = precast_k10, fill = precast_k10)) + geom_density_ridges()+ theme_ridges()
