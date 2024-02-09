library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)
library(nnSVG)
library(cowplot)

spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered_hdf5"
)
out_dir <- here("processed-data", "05_harmony_BayesSpace")
plot_dir <- here("plots", "05_harmony_BayesSpace")
nn_methods = c("default", "precast")

hvg_scran_path <- here(
    "processed-data", "05_harmony_BayesSpace", "top.hvgs.Rdata"
)
sig_cutoff <- 0.05
best_sample_id <- "Br8492"

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)

top_svgs = list()

for (method in c("default", "precast")) {
    if (method == "default") {
        nn_out_dir = file.path(out_dir, 'nnSVG_out')
        nn_plot_dir = file.path(plot_dir, 'nnSVG')
    } else {
        nn_out_dir = file.path(out_dir, 'nnSVG_precast_out')
        nn_plot_dir = file.path(plot_dir, 'nnSVG_precast')
    }

    ############################################################################
    #   Gather nnSVG results across samples and compute summary metrics
    ############################################################################

    nn_out_list <- list()
    for (sample_id in unique(spe$sample_id)) {
        nn_out_list[[sample_id]] <- file.path(
                nn_out_dir, paste0(sample_id, ".csv")
            ) |>
            read.csv() |>
            as_tibble() |>
            mutate(sample_id = sample_id)
    }
    nn_out <- do.call(rbind, nn_out_list)

    #   Compute summary metrics across donors: proportion of samples where the
    #   gene was significant; proportion in the top 100 ranks; average rank
    num_samples <- length(unique(spe$sample_id))
    nn_out_summary <- nn_out |>
        group_by(gene_id) |>
        summarize(
            nnsvg_prop_sig_adj = sum(padj < sig_cutoff) / num_samples,
            nnsvg_prop_top_100 = sum(rank < 100) / num_samples,
            nnsvg_avg_rank = mean(rank)
        ) |>
        #   Compute a rank of average ranks
        mutate(
            nnsvg_avg_rank_rank = match(nnsvg_avg_rank, sort(nnsvg_avg_rank))
        )

    write_csv(nn_out_summary, file.path(nn_out_dir, "summary_across_samples.csv"))

    ############################################################################
    #   Visually examine top-ranked SVGs
    ############################################################################

    #   Order rank of genes of those where all samples had statistically
    #   significant spatial variability
    top_svgs[[method]] <- nn_out_summary |>
        filter(nnsvg_prop_sig_adj == 1) |>
        arrange(nnsvg_avg_rank_rank) |>
        mutate(
            symbol = rowData(spe)[
                match(gene_id, rowData(spe)$gene_id), "gene_name"
            ],
            method_name = method
        )
    
    #   Export top 100 SVGs for this method
    top_svgs[[method]] |>
        select(gene_id, symbol) |>
        slice_head(n = 100) |>
        write_csv(file.path(nn_out_dir, "top_100_SVGs.csv"))

    #   Plot the top 50 genes for one sample
    num_genes <- 50
    plot_list <- list()
    for (i in 1:num_genes) {
        plot_list[[i]] <- spot_plot(
            spe,
            sample_id = best_sample_id,
            title = paste(
                best_sample_id, top_svgs[[method]]$symbol[i], sep = "_"
            ),
            var_name = top_svgs[[method]]$gene_id[i],
            is_discrete = FALSE,
            minCount = 0
        )
    }

    pdf(file.path(nn_plot_dir, sprintf("nnSVG_top_SVGs_%s.pdf", best_sample_id)))
    print(plot_list)
    dev.off()

    #   Plot a grid of several samples and several top-ranked genes
    num_genes <- 5
    num_samples <- 5
    plot_list <- list()
    counter <- 1
    for (i in 1:num_genes) {
        for (j in 1:num_samples) {
            this_sample_id <- unique(spe$sample_id)[j]
            plot_list[[counter]] <- spot_plot(
                spe,
                sample_id = this_sample_id,
                title = paste(
                    this_sample_id, top_svgs[[method]]$symbol[i], sep = "_"
                ),
                var_name = top_svgs[[method]]$gene_id[i],
                is_discrete = FALSE,
                minCount = 0
            )
            counter <- counter + 1
        }
    }

    pdf(
        file.path(nn_plot_dir, "nnSVG_grid.pdf"),
        width = 7 * num_samples,
        height = 7 * num_genes
    )
    print(plot_grid(plotlist = plot_list, ncol = num_samples))
    dev.off()
}

################################################################################
#   Form a tibble of top HVGs (from two methods: scran::getTopHVGs and based on
#   binomial deviance)
################################################################################

top_hvgs = list()

#   Load HVGs from scran::getTopHVGs, add gene symbols, and take just
#   significant (after adjustment) genes
load(hvg_scran_path, verbose = TRUE)
top_hvgs[['scran']] <- tibble(
    gene_id = top.hvgs.fdr5,
    gene_name = rowData(spe)[
        match(gene_id, rowData(spe)$gene_id), "gene_name"
    ],
    method_name = 'scran'
)

if (any(is.na(top_hvgs$gene_name))) {
    stop("Some HVGs not in SpatialExperiment")
}

#   Prepare highly deviant genes
top_hvgs[['deviance']] = rowData(spe) |>
    as_tibble() |>
    arrange(desc(binomial_deviance)) |>
    select(gene_id, gene_name) |>
    mutate(method_name = 'deviance')

#   Drop mitochondrial genes (as was done for SVGs, considering that expression
#   variation in these genes represents technical noise)
top_hvgs = do.call(rbind, top_hvgs) |>
    filter(!str_detect(gene_name, "^MT-"))


#   Also just form one tibble of top SVGs, like that for the HVGs
top_svgs = do.call(rbind, top_svgs)

################################################################################
#   Compare each combination of SVGs and HVGs (2 methods for each)
################################################################################

num_points <- 100
overlap_df <- tibble(
    #   Sample [num_points] different numbers of genes, linearly spaced between
    #   0 and as many top genes are shared
    num_genes = as.integer(
        1:num_points * min(nrow(top_svgs), nrow(top_hvgs)) / num_points
    ),
    prop_overlap = NA
)

#   At each number of genes, compute the proportion of SVGs that are also HVGs
#   at the same cutoff
for (i in 1:num_points) {
    overlap_df[i, "prop_overlap"] <- mean(
        head(top_svgs$gene_id, overlap_df$num_genes[i]) %in%
            head(top_hvgs$gene_id, overlap_df$num_genes[i])
    )
}

#   Plot how this proportion changes with number of genes
p <- ggplot(overlap_df) +
    geom_line(aes(x = num_genes, y = prop_overlap)) +
    scale_y_continuous(limits = c(0, max(overlap_df$prop_overlap)))
pdf(file.path(plot_dir, "HVG_SVG_prop_overlap.pdf"))
print(p)
dev.off()

#   Export top 100 HVGs for this method
top_hvgs |>
    dplyr::rename(symbol = gene_name) |>
    slice_head(n = 100) |>
    write_csv(file.path(dirname(hvg_path), "top_100_HVGs.csv"))

session_info()
