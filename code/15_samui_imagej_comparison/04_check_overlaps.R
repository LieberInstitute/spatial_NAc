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
final_steps = c('imagej', 'samui')
k_values = c(2, 10, 19, 28)

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

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$capture_area = spe$sample_id_original

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
    mutate(
        array_row = ifelse(
            final_step == 'imagej', array_row_imagej, array_row_samui
        ),
        array_col = ifelse(
            final_step == 'imagej', array_col_imagej, array_col_samui
        ),
        exclude_overlapping = ifelse(
            final_step == 'imagej', exclude_overlapping_imagej, exclude_overlapping_samui
        )
    ) |>
    select(-matches('^array_(row|col)_(imagej|samui)$'))

stopifnot(!any(is.na(cluster_df)))
