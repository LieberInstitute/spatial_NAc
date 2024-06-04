library(getopt)
library(sessioninfo)
library(here)
library(SpatialExperiment)
library(tidyverse)
library(HDF5Array)
library(visiumStitched)

spe_dir = here('processed-data', '15_samui_imagej_comparison', 'spe')
bs_cluster_path = here(
    'processed-data', '15_samui_imagej_comparison', 'bayesspace_out',
    'k%s_%s.csv'
)
precast_cluster_path = here(
    'processed-data', '15_samui_imagej_comparison', 'precast_out',
    'PRECAST_k%s_%s.csv'
)
plot_dir = here('plots', '15_samui_imagej_comparison')
final_steps = c('imagej', 'samui')
best_sample_id = 'Br8492'
k_values = c(2, 10, 19, 28)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Read in cluster assignments
################################################################################

message('Reading in cluster assignments from PRECAST and BayesSpace...')
cluster_df_list = list()
i = 1

#   Individually read in PRECAST and BayesSpace results for all values of k,
#   before and after refinement with Samui
for (this_k in k_values) {
    for (this_final_step in final_steps) {
        precast_df = read_csv(
                sprintf(precast_cluster_path, this_k, this_final_step),
                show_col_types = FALSE
            ) |>
            dplyr::rename(cluster_assignment = cluster) |>
            dplyr::select(key, cluster_assignment) |>
            mutate(
                final_step = factor(this_final_step, levels = final_steps),
                k = factor(this_k, levels = k_values),
                method = factor('PRECAST', levels = c('PRECAST', 'BayesSpace'))
            )
        bs_df = read_csv(
                sprintf(bs_cluster_path, this_k, this_final_step),
                show_col_types = FALSE
            ) |>
            mutate(
                final_step = factor(this_final_step, levels = final_steps),
                k = factor(this_k, levels = k_values),
                method = factor(
                    'BayesSpace', levels = c('PRECAST', 'BayesSpace')
                )
            )
        cluster_df_list[[i]] = rbind(precast_df, bs_df)
        i = i + 1
    }
}

#   Combine into a single tibble
cluster_df = do.call(rbind, cluster_df_list)

################################################################################
#   Read in SPE, add overlap info, and merge with cluster assignments
################################################################################

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$capture_area = spe$sample_id_original

#   Add overlap info for both forms of array coordinates. Really should've done
#   this in 01_build_spe.R...
for (final_step in final_steps) {
    message(sprintf('Adding overlap info for final_step=%s...', final_step))

    #   Find overlap info with the array coordinates given by 'final_step'
    spe$array_row_transformed = spe[[paste0('array_row_', final_step)]]
    spe$array_col_transformed = spe[[paste0('array_col_', final_step)]]
    spe = add_overlap_info(spe, 'sum_umi')

    #   Rename the columns with 'final_step' that were added by
    #   'add_overlap_info()'
    spe[[paste0('exclude_overlapping_', final_step)]] = spe$exclude_overlapping
    spe[[paste0('overlap_key_', final_step)]] = spe$overlap_key

    #   Remove the extra unnecessary columns added when finding overlap info
    temp = colData(spe) |>
        as_tibble() |>
        select(
            -c(
                matches('^array_(row|col)_transformed'), exclude_overlapping,
                overlap_key
            )
        ) |>
        DataFrame()
    rownames(temp) = colnames(spe)
    colData(spe) = temp
}

col_data = colData(spe) |>
    as_tibble() |>
    select(
        matches('^array_(row|col)_(imagej|samui)$'),
        matches('^(exclude_overlapping|overlap_key)_(imagej|samui)$'),
        capture_area, key, donor
    )

cluster_df = cluster_df |>
    left_join(col_data, by = "key") |>
    #   Only take the array coordinates and overlap info related to the final
    #   step in each row
    mutate(
        array_row = ifelse(
            final_step == 'imagej', array_row_imagej, array_row_samui
        ),
        array_col = ifelse(
            final_step == 'imagej', array_col_imagej, array_col_samui
        ),
        exclude_overlapping = ifelse(
            final_step == 'imagej', exclude_overlapping_imagej, exclude_overlapping_samui
        ),
        overlap_key = ifelse(
            final_step == 'imagej', overlap_key_imagej, overlap_key_samui
        )
    ) |>
    select(
        -c(
            matches('^array_(row|col)_(imagej|samui)$'),
            matches('^(exclude_overlapping|overlap_key)_(imagej|samui)$')
        )
    )

################################################################################
#   Plot clustering results spatially
################################################################################

for (this_final_step in final_steps) {
    for (this_k in k_values) {
        for (this_method in c('PRECAST', 'BayesSpace')) {
            #   Add clustering results to the SPE for this combination of
            #   parameters
            temp = colData(spe) |>
                as_tibble() |>
                select(key, sample_id) |>
                left_join(
                    cluster_df |>
                        filter(
                            final_step == this_final_step,
                            k == this_k,
                            method == this_method
                        ) |>
                        select(
                            key, cluster_assignment, exclude_overlapping,
                            array_row, array_col
                        ),
                    by = "key"
                ) |>
                DataFrame()
            rownames(temp) = colnames(spe)
            colData(spe) = temp

            p = spot_plot(
                spe, sample_id = best_sample_id,
                var_name = 'cluster_assignment', is_discrete = TRUE
            )
            pdf(
                file.path(
                    plot_dir,
                    sprintf(
                        '%s_%s_%s_%s.pdf',
                        best_sample_id, this_method, this_k, this_final_step
                    )
                )
            )
            print(p)
            dev.off()
        }
    }
}

################################################################################
#   Merge and explore data at overlaps
################################################################################

overlap_df = cluster_df |>
    #   Only take spots overlapping other spots   
    filter(overlap_key != "") |>
    #   Add cluster assignment of the overlapping spot
    group_by(method) |>
    mutate(
        overlap_cluster_assignment = cluster_assignment[match(overlap_key, key)]
    ) |>
    ungroup() |>
    #   Remove duplicate pairs and uncommon cases where the overlapping spot
    #   was dropped before clustering
    filter(!(key %in% overlap_key), !is.na(overlap_cluster_assignment))

stopifnot(!any(is.na(overlap_df)))

p = overlap_df |>
    group_by(final_step, k, method, donor) |>
    summarize(
        agreement = mean(cluster_assignment == overlap_cluster_assignment)
    ) |>
    ungroup() |>
    ggplot(mapping = aes(x = k, y = agreement, color = final_step)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) +
        facet_wrap(~method) +
        labs(
            y = "Proportion agreement at overlaps",
            color = "Final alignment step"
        ) +
        theme_bw(base_size = 18)

pdf(file.path(plot_dir, 'all_overlaps.pdf'), width = 10, height = 7)
print(p)
dev.off()

session_info()
