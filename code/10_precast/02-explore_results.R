library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)

out_path = here('processed-data', '10_precast', 'PRECAST_k%s.csv')
spe_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
wm_genes = c('MBP', 'GFAP', 'PLP1', 'AQP4')

spe = loadHDF5SummarizedExperiment(spe_dir)

################################################################################
#   Import PRECAST results and append to colData
################################################################################

result_list = list()
for (K in 2:28) {
    result_list[[K]] = sprintf(out_path, K) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster)) |>
        mutate(k = K)
}
precast_results = do.call(rbind, result_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "precast_k"
    )

temp = colnames(spe)
colData(spe) = colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) = temp

################################################################################
#   Check that k = 2 strongly correlates with WM vs GM boundary
################################################################################

#   For each spot, average expression Z-scores across all white matter genes.
#   The idea is that the combined expression of several white matter genes
#   should form a more precise marker of white matter than one gene
stopifnot(all(wm_genes %in% rowData(spe)$gene_name))
a <- assays(spe)$logcounts[match(wm_genes, rowData(spe)$gene_name),]
gene_z <- (a - rowMeans(a)) / rowSds(a)
spe$z_score <- colMeans(gene_z, na.rm = TRUE)

#   The averaged Z-score of expression across white matter genes should strongly
#   correlate with the cluster assignment at k=2
print("Overall correlation of WM gene expression with k=2 cluster assignment:")
print(cor(spe$z_score, spe$precast_k2, use = "complete.obs"))

print("Correlation broken down by donor:")
colData(spe) |>
    as_tibble() |>
    group_by(sample_id) |>
    summarize(cor_value = cor(z_score, precast_k2, use = "complete.obs")) |>
    print()

#   Note that correlations are consistently large and positive, indicating that
#   cluster 1 roughly corresponds to gray matter, and 2 roughly to white matter
