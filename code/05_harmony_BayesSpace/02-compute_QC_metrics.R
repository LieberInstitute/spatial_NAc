library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(scran)
library(scater)
library(scry)
library(HDF5Array)
library(bluster)
library(SpotSweeper)
library(pheatmap)
library(ggspavis)
library(ggplot2)
library(ggExtra)

# Read in the raw SPE object
code_dir = here("code", "05_harmony_BayesSpace")
source(paste0(code_dir, "/localOutlier_2.R"))

processed_dir = here("processed-data", "05_harmony_BayesSpace")
raw_in_path = here('processed-data', '05_harmony_BayesSpace', '01-build_spe', 'spe_raw.rds')
QC_added_ordinary_path = here(
    'processed-data', '05_harmony_BayesSpace', '02-compute_QC_metrics','spe_with_QC_metrics.rds'
)
QC_added_hdf5_dir = here(
    'processed-data', '05_harmony_BayesSpace', '02-compute_QC_metrics', 'spe_with_QC_metrics_hdf5'
)
plot_dir = here('plots', '05_harmony_BayesSpace', '02-compute_QC_metrics')

num_cores = Sys.getenv('SLURM_CPUS_ON_NODE')
set.seed(0)

spe = readRDS(raw_in_path)
cat("Initial number of spots:", dim(spe)[2], "\n")

# Add the number of estimated cells per spot using Vistoseg
sample_info_path = here('processed-data', 'VistoSeg', 'VistoSeg_inputs.csv')
sample_info = read.csv(sample_info_path) |> as_tibble()

counts_list = list()
for (i in 1:nrow(sample_info)) {
    counts_list[[i]] = read.csv(
                file.path(sample_info$spaceranger_dir[i],
                'tissue_spot_counts.csv')
            ) |>
            as_tibble() |>
            mutate(sample_id = sample_info$sample_id[i])
}

counts = do.call(rbind, counts_list) |>
    mutate(key = paste(barcode, sample_id, sep = '_')) |>
    filter(key %in% spe$key)
counts <- data.frame(counts)
counts <- counts[match(spe$key, counts$key), ]
spe$Nmask_dark_blue <- counts$Nmask_dark_blue
spe$Pmask_dark_blue <- counts$Pmask_dark_blue
spe$CNmask_dark_blue <- counts$CNmask_dark_blue

df <- colData(spe)
df <- data.frame(df)
ggplot(df, aes(x = sum_umi, fill = in_tissue)) + geom_density(alpha=0.25) + scale_x_continuous(trans='log10') + theme_classic()

# Add manual annotations by Svitlana
manual_anno_path = here('raw-data', 'manual_annotations')
unique_samples <- unique(spe$sample_id)
manual_annotations <- lapply(unique_samples, function(sample_id){
    cat(sample_id, "\n")
    df <- read.csv(file.path(manual_anno_path, paste0(sample_id, ".csv")))
    colnames(df) <- c("sample", "spot", "annotation")
    df
})
manual_annotations <- do.call(rbind, manual_annotations)



spe <- spe[, (colSums(assays(spe)$counts) > 0) & spe$in_tissue]
cat("Number of spots after preliminary QC:", dim(spe)[2], "\n")

## Metrics QC
metrics_qc <- function(spe) {

    qc_df <- data.frame(
        log2sum = log2(spe$sum_umi),
        log2detected = log2(spe$sum_gene),
        subsets_Mito_percent = spe$expr_chrM_ratio*100,
        sample_id = spe$sample_id_original
    )

    qcfilter <- DataFrame(
        low_lib_size = isOutlier(qc_df$log2sum, type = "lower", log = TRUE, batch = qc_df$sample_id),
        low_n_features = isOutlier(qc_df$log2detected, type = "lower", log = TRUE, batch = qc_df$sample_id),
        high_subsets_Mito_percent = isOutlier(qc_df$subsets_Mito_percent, type = "higher", batch = qc_df$sample_id)
    )
    qcfilter$discard <- (qcfilter$low_lib_size | qcfilter$low_n_features) | qcfilter$high_subsets_Mito_percent


    spe$scran_low_lib_size_low_mito <- factor(qcfilter$low_lib_size & qc_df$subsets_Mito_percent < 0.5, levels = c("TRUE", "FALSE"))


    spe$scran_discard <-
        factor(qcfilter$discard, levels = c("TRUE", "FALSE"))
    spe$scran_low_lib_size <-
        factor(qcfilter$low_lib_size, levels = c("TRUE", "FALSE"))
    spe$scran_low_n_features <-
        factor(qcfilter$low_n_features, levels = c("TRUE", "FALSE"))
    spe$scran_high_subsets_Mito_percent <-
        factor(qcfilter$high_subsets_Mito_percent, levels = c("TRUE", "FALSE"))

    ## Find edge spots
    spots <- data.frame(
        row = spe$array_row,
        col = spe$array_col,
        sample_id = spe$sample_id_original
    )

    edge_spots_row <- group_by(spots, sample_id, row) %>% summarize(min_col = min(col), max_col = max(col))
    edge_spots_col <- group_by(spots, sample_id, col) %>% summarize(min_row = min(row), max_row = max(row))

    spots <- left_join(spots, edge_spots_row) %>% left_join(edge_spots_col)
    spots$edge_spots <- with(spots, row == min_row | row == max_row | col == min_col | col == max_col)

    spots$row_distance <- with(spots, pmin(abs(row - min_row), abs(row - max_row)))
    spots$col_distance <- with(spots, pmin(abs(col - min_col), abs(col - max_col)))
    ## spots$edge_distance <- with(spots, sqrt(row_distance^2 + col_distance^2))
    ## The above is from:
    ## sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2)
    ## but it was wrong, here's a case the the smallest distance is on the column:
    ## sqrt(0^2 + col_distance^2) = col_distance
    spots$edge_distance <- with(spots, pmin(row_distance, col_distance))


    spe$edge_spots <- factor(spots$edge_spots, levels = c("TRUE", "FALSE"))
    spe$edge_distance <- spots$edge_distance


    spe$scran_low_lib_size_edge <- factor(qcfilter$low_lib_size & spots$edge_distance < 1, levels = c("TRUE", "FALSE"))

    return(spe)
}

spe <- metrics_qc(spe)

# Check which spots are found to be outliers based on SpotSweeper
spe$array_row <- spe$array_row_original
spe$array_col <- spe$array_col_original
spe$sample_id_original <- as.character(spe$sample_id_original)
# Find SpotSweeper outliers
feats <- c("sum_umi", "sum_gene", "expr_chrM_ratio")
spe <- localOutliers_2(spe, features = feats, n_neighbors=18, 
                    data_output=TRUE,
                    method="multivariate", samples = "sample_id_original", n_cores = num_cores)

spe$array_row <- spe$array_row_transformed
spe$array_col <- spe$array_col_transformed
spe$sample_id_original <- as.factor(spe$sample_id_original)


# Save data with global and local outlier data
message(Sys.time(), " - Saving HDF5-backed filtered spe")
spe = saveHDF5SummarizedExperiment(
    spe, dir = QC_added_hdf5_dir, replace = TRUE
)

message(Sys.time(), " - Saving ordinary filtered spe")
spe = realize(spe)
saveRDS(spe, QC_added_ordinary_path)

pdf(width = 8, height = 8, paste0(plot_dir, "/sum_umi_spot_plot.pdf"))
spot_plot(spe, "Br2720", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br2743", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br3942", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6423", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6432", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6471", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6522", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8325", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8492", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8667", var_name = "sum_umi_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
dev.off()

pdf(width = 8, height = 8, paste0(plot_dir, "/sum_gene_spot_plot.pdf"))
spot_plot(spe, "Br2720", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br2743", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br3942", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6423", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6432", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6471", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br6522", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8325", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8492", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
spot_plot(spe, "Br8667", var_name = "sum_gene_log2", is_discrete = FALSE, spatial = TRUE, assayname = "counts")
dev.off()

pdf(width = 8, height = 8, paste0(plot_dir, "/mito_ratio_spot_plot.pdf"))
spot_plot(spe, "Br2720", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br2743", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br3942", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br6423", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br6432", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br6471", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br6522", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br8325", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br8492", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
spot_plot(spe, "Br8667", var_name = "expr_chrM_ratio", is_discrete = FALSE, spatial = TRUE, assayname = "counts", minCount = 0)
dev.off()

session_info()

