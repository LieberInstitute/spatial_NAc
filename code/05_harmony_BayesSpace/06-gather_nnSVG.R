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
hvg_path = here("processed-data", "05_harmony_BayesSpace", "top.hvgs.Rdata")
plot_dir = here('plots', '05_harmony_BayesSpace')
sig_cutoff = 0.05
best_sample_id = 'Br8492'

spe = loadHDF5SummarizedExperiment(spe_dir)

################################################################################
#   Gather nnSVG results across samples and compute summary metrics
################################################################################

nn_out_list = list()
for (sample_id in unique(spe$sample_id)) {
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

################################################################################
#   Visually examine top-ranked SVGs
################################################################################

top_genes = nn_out_summary |>
    filter(nnsvg_avg_rank_rank <= 5) |>
    pull(gene_id)

plot_list = lapply(
    top_genes,
    function(gene) {
        spot_plot(
            spe,
            sample_id = best_sample_id,
            var_name = gene,
            is_discrete = FALSE,
            minCount = 0
        )
    }
)

pdf(file.path(plot_dir, 'nnSVG_top_SVGs.pdf'))
print(plot_list)
dev.off()

################################################################################
#   Update SpatialExperiment with nnSVG results
################################################################################

#   Add summary metrics from nnSVG to the SpatialExperiment
rowData(spe) = rowData(spe) |>
    as_tibble() |>
    left_join(nn_out_summary, by = 'gene_id') |>
    column_to_rownames('gene_id') |>
    DataFrame()
quickResaveHDF5SummarizedExperiment(spe)

nn_out_summary |>
    filter(nnsvg_avg_rank_rank <= 0.1 * nrow(spe))

load(hvg_path)
session_info()
