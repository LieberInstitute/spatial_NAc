library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)
library(nnSVG)

spe_dir = here(
    'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)
nn_out_dir = here('processed-data', '05_harmony_BayesSpace', 'nnSVG_out')
sig_cutoff = 0.05

spe = loadHDF5SummarizedExperiment(spe_dir)

nn_out_list = list()
for (sample_id in unique(spe$sample_id)[2:10]) {
    nn_out_list[[sample_id]] = file.path(nn_out_dir, paste0(sample_id, '.csv')) |>
        read.csv() |>
        as_tibble()
}
nn_out = do.call(rbind, nn_out_list)

#   Compute summary metrics across donors: proportion of samples where the gene
#   was significant; proportion in the top 100 ranks; average rank
nn_out_summary = nn_out |>
    group_by(gene_id) |>
    summarize(
        nnsvg_prop_sig_adj = mean(padj < sig_cutoff),
        nnsvg_prop_top_100 = mean(rank < 100),
        nnsvg_avg_rank = mean(rank)
    ) |>
    #   Compute a rank of average ranks
    mutate(nnsvg_avg_rank_rank = match(nnsvg_avg_rank, sort(nnsvg_avg_rank)))

#   Add summary metrics from nnSVG to the SpatialExperiment
rowData(spe) = rowData(spe) |>
    as_tibble() |>
    left_join(nn_out_summary, by = 'gene_id') |>
    column_to_rownames('gene_id') |>
    DataFrame()
quickResaveHDF5SummarizedExperiment(spe)

session_info()
