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

################################################################################
#   For each method of running nnSVG, gather results across samples and
#   compute summary metrics
################################################################################

top_svgs = list()
for (method in c("default", "precast")) {
    if (method == "default") {
        nn_out_dir = file.path(out_dir, 'nnSVG_out')
    } else {
        nn_out_dir = file.path(out_dir, 'nnSVG_precast_out')
    }

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

    #   Export the raw summary, not dropping any results according to rank
    #   or statistical significance
    write_csv(nn_out_summary, file.path(nn_out_dir, "summary_across_samples.csv"))

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
}

################################################################################
#   Form a tibble of top HVGs (from two methods: scran::getTopHVGs and based on
#   binomial deviance)
################################################################################

top_hvgs = list()

#   Load significant (after adjustment) HVGs from scran::getTopHVGs
load(hvg_scran_path, verbose = TRUE)
top_hvgs[['scran']] <- tibble(
    gene_id = top.hvgs.fdr5,
    method_name = 'scran'
) |>
    mutate(rank = row_number())

#   Prepare highly deviant genes
top_hvgs[['deviance']] = rowData(spe) |>
    as_tibble() |>
    arrange(desc(binomial_deviance)) |>
    select(gene_id) |>
    mutate(
        method_name = 'deviance',
        rank = row_number()
    )

top_hvgs = do.call(rbind, top_hvgs) |>
    mutate(variation_type = 'HVG')

################################################################################
#   Compute both methods of HVGs and SVGs into one tibble, 'top_genes'
################################################################################

#   Also just form one tibble of top SVGs, like that for the HVGs
top_svgs = do.call(rbind, top_svgs) |>
    group_by(method_name) |>
    arrange(nnsvg_avg_rank_rank) |>
    mutate(
        variation_type = 'SVG',
        #   Re-rank, since not-significant SVGs were dropped earlier
        rank = row_number()
    ) |>
    ungroup() |>
    select(colnames(top_hvgs))
    
top_genes = rbind(top_hvgs, top_svgs) |>
    #   Add gene symbol
    mutate(
        gene_name = rowData(spe)[
            match(gene_id, rowData(spe)$gene_id), "gene_name"
        ]
    ) |>
    #   Drop mitochondrial genes (already done for SVGs), considering that
    #   expression variation in these genes represents technical noise
    filter(!str_detect(gene_name, "^MT-")) |>
    #   Re-rank genes after dropping
    group_by(method_name, variation_type) |>
    arrange(rank) |>
    mutate(rank = row_number()) |>
    ungroup()

if (any(is.na(top_genes$gene_name))) {
    stop("Some top genes not in SpatialExperiment")
}

write_csv(top_genes, file.path(out_dir, 'top_variable_genes.csv'))

################################################################################
#   Compare each combination of SVGs and HVGs (2 methods for each)
################################################################################

svg_methods = top_genes |>
    filter(variation_type == 'SVG') |>
    pull(method_name) |>
    unique()
hvg_methods = top_genes |>
    filter(variation_type == 'HVG') |>
    pull(method_name) |>
    unique()

#   Get the minimum number of top genes from any method
most_shared = top_genes |>
    group_by(method_name, variation_type) |>
    summarize(total_genes = max(rank)) |>
    pull(total_genes) |>
    min()

overlap_df_list = list()
for (svg_method in svg_methods) {
    for (hvg_method in hvg_methods) {
        these_svgs = top_genes |>
            filter(method_name == svg_method, variation_type == 'SVG')
        these_hvgs = top_genes |>
            filter(method_name == hvg_method, variation_type == 'HVG')
        combined_method = paste(svg_method, hvg_method, sep = '_')

        num_points <- 100
        overlap_df_list[[combined_method]] <- tibble(
            #   Sample [num_points] different numbers of genes, linearly spaced
            #   between 0 and as many top genes are shared
            num_genes = as.integer(
                1:num_points * most_shared / num_points
            ),
            prop_overlap = NA,
            method_name = combined_method
        )

        #   At each number of genes, compute the proportion of SVGs that are
        #   also HVGs at the same cutoff
        for (i in 1:num_points) {
            overlap_df_list[[combined_method]][i, "prop_overlap"] <- mean(
                head(
                    these_svgs$gene_id,
                    overlap_df_list[[combined_method]]$num_genes[i]
                ) %in%
                head(
                    these_hvgs$gene_id,
                    overlap_df_list[[combined_method]]$num_genes[i]
                )
            )
        }
    }
}

overlap_df = do.call(rbind, overlap_df_list)

#   Plot how this proportion changes with number of genes
p <- ggplot(overlap_df) +
    geom_line(aes(x = num_genes, y = prop_overlap)) +
    facet_wrap(~method_name) +
    scale_y_continuous(limits = c(0, max(overlap_df$prop_overlap))) +
    labs(x = 'Top N genes', y = 'Prop. Shared Genes')
pdf(file.path(plot_dir, "HVG_SVG_prop_overlap.pdf"))
print(p)
dev.off()

################################################################################
#   Plot and export top variable genes for each method of HVGs and SVGs
################################################################################

for (this_method in unique(top_genes$method_name)) {
    these_genes = top_genes |>
        filter(method_name == this_method)
    
    method_title = sprintf(
        '%s_%ss', this_method, unique(these_genes$variation_type)
    )

    #   Export top 100 genes for this method
    these_genes |>
        select(gene_id, gene_name, rank) |>
        dplyr::rename(symbol = gene_name) |>
        slice_head(n = 100) |>
        write_csv(file.path(out_dir, sprintf("top_100_%s.csv", method_title)))
    
    #   Plot the top 50 genes for one sample
    num_genes <- 50
    plot_list <- list()
    for (i in 1:num_genes) {
        plot_list[[i]] <- spot_plot(
            spe,
            sample_id = best_sample_id,
            title = paste(
                best_sample_id, these_genes$gene_name[i], sep = "_"
            ),
            var_name = these_genes$gene_id[i],
            is_discrete = FALSE,
            minCount = 0
        )
    }

    pdf(
        file.path(
            plot_dir,
            sprintf("top_50_%s_%s.pdf", method_title, best_sample_id)
        )
    )
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
                    this_sample_id, these_genes$gene_name[i], sep = "_"
                ),
                var_name = these_genes$gene_id[i],
                is_discrete = FALSE,
                minCount = 0
            )
            counter <- counter + 1
        }
    }

    pdf(
        file.path(plot_dir, paste(method_title, "grid.pdf", sep = '_')),
        width = 7 * num_samples,
        height = 7 * num_genes
    )
    print(plot_grid(plotlist = plot_list, ncol = num_samples))
    dev.off()
}

session_info()
