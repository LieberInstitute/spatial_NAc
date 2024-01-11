library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)
library(nnSVG)
library(cowplot)

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
        as_tibble() |>
        mutate(sample_id = sample_id)
}
nn_out = do.call(rbind, nn_out_list)

#   Compute summary metrics across donors: proportion of samples where the gene
#   was significant; proportion in the top 100 ranks; average rank
num_samples = length(unique(spe$sample_id))
nn_out_summary = nn_out |>
    group_by(gene_id) |>
    summarize(
        nnsvg_prop_sig_adj = sum(padj < sig_cutoff) / num_samples,
        nnsvg_prop_top_100 = sum(rank < 100) / num_samples,
        nnsvg_avg_rank = mean(rank)
    ) |>
    #   Compute a rank of average ranks
    mutate(nnsvg_avg_rank_rank = match(nnsvg_avg_rank, sort(nnsvg_avg_rank)))

write_csv(nn_out_summary, file.path(nn_out_dir, 'summary_across_samples.csv'))

################################################################################
#   Visually examine top-ranked SVGs
################################################################################

#   Order rank of genes of those where all samples had statistically
#   significant spatial variability
top_genes = nn_out_summary |>
    filter(nnsvg_prop_sig_adj == 1) |>
    arrange(nnsvg_avg_rank_rank) |>
    mutate(
        symbol = rowData(spe)[match(gene_id, rowData(spe)$gene_id), 'gene_name']
    )

#   Plot the top 50 genes for one sample
num_genes = 50
plot_list = list()
for (i in 1:num_genes) {
    plot_list[[i]] = spot_plot(
        spe,
        sample_id = best_sample_id,
        title = paste(best_sample_id, top_genes$symbol[i], sep = '_'),
        var_name = top_genes$gene_id[i],
        is_discrete = FALSE,
        minCount = 0
    )
}

pdf(file.path(plot_dir, sprintf('nnSVG_top_SVGs_%s.pdf', best_sample_id)))
print(plot_list)
dev.off()

#   Plot a grid of several samples and several top-ranked genes
num_genes = 5
num_samples = 5
plot_list = list()
counter = 1
for (i in 1:num_genes) {
    for (j in 1:num_samples) {
        this_sample_id = unique(spe$sample_id)[j]
        plot_list[[counter]] = spot_plot(
            spe,
            sample_id = this_sample_id,
            title = paste(this_sample_id, top_genes$symbol[i], sep = '_'),
            var_name = top_genes$gene_id[i],
            is_discrete = FALSE,
            minCount = 0
        )
        counter = counter + 1
    }
}

pdf(
    file.path(plot_dir, 'nnSVG_grid.pdf'),
    width = 7 * num_samples,
    height = 7 * num_genes
)
print(plot_grid(plotlist = plot_list, ncol = num_samples))
dev.off()

################################################################################
#   Compare with HVGs (as ranked by 'getTopHVGs')
################################################################################

#   Load HVGs, add gene symbols, and take just significant (after adjustment),
#   non-mitochondrial genes
load(hvg_path, verbose = TRUE)
top_genes_hvg = tibble(
        gene_id = top.hvgs.fdr5,
        gene_name = rowData(spe)[
            match(gene_id, rowData(spe)$gene_id), 'gene_name'
        ]
    ) |>
    filter(!str_detect(gene_name, '^MT-'))

if (any(is.na(top_genes_hvg$gene_name))) {
    stop("Some HVGs not in SpatialExperiment")
}

num_points = 100
overlap_df = tibble(
    #   Sample [num_points] different numbers of genes, linearly spaced between
    #   0 and as many top genes are shared
    num_genes = as.integer(
        1:num_points * min(nrow(top_genes), nrow(top_genes_hvg)) / num_points
    ),
    prop_overlap = NA
)

#   At each number of genes, compute the proportion of SVGs that are also HVGs
#   at the same cutoff 
for (i in 1:num_points) {
    overlap_df[i, 'prop_overlap'] = mean(
        head(top_genes$gene_id, overlap_df$num_genes[i]) %in%
        head(top_genes_hvg$gene_id, overlap_df$num_genes[i])
    )
}

#   Plot how this proportion changes with number of genes
p = ggplot(overlap_df) +
    geom_line(aes(x = num_genes, y = prop_overlap)) +
    scale_y_continuous(limits = c(0, max(overlap_df$prop_overlap)))
pdf(file.path(plot_dir, 'HVG_SVG_prop_overlap.pdf'))
print(p)
dev.off()

################################################################################
#   Export CSV of top SVGs and HVGs
################################################################################

top_genes |>
    select(gene_id, symbol) |>
    slice_head(n = 100) |>
    write_csv(file.path(dirname(hvg_path), 'top_100_SVGs.csv'))

top_genes_hvg |>
    dplyr::rename(symbol = gene_name) |>
    slice_head(n = 100) |>
    write_csv(file.path(dirname(hvg_path), 'top_100_HVGs.csv'))

session_info()
