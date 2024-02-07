library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(purrr)

out_path <- here("processed-data", "10_precast", "PRECAST_k%s.csv")
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "spe_filtered_hdf5"
)
plot_dir <- here("plots", "10_precast")

#   List of genes provided by Svitlana, named by the subregion they're markers
#   for
subregion_genes <- list(
    "shell" = c(
        "TAC3", "TAC1", "CALB2", "GRIA4", "CPNE4", "GREB1L", "ARHGAP6", "OPRK1"
    ),
    "core" = toupper(
        c(
            "CALB1", "HPCAL4", "ZDBF2", "LAMP5", "Scn4b", "Rgs9", "Ppp3ca",
            "Oprd1", "Adora2a", "Pde1a", "Peg10", "Dlk1", "Htr2c"
        )
    ),
    "white_matter" = c("MBP", "GFAP", "PLP1", "AQP4"),
    "D1_islands" = c("OPRM1", "CHST9")
)

dir.create(plot_dir, showWarnings = FALSE)

spe <- loadHDF5SummarizedExperiment(spe_dir)
stopifnot(all(unlist(subregion_genes) %in% rowData(spe)$gene_name))

################################################################################
#   Import PRECAST results and append to colData
################################################################################

result_list <- list()
for (K in 2:28) {
    result_list[[K]] <- sprintf(out_path, K) |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster)) |>
        mutate(k = K)
}
precast_results <- do.call(rbind, result_list) |>
    pivot_wider(
        values_from = cluster, names_from = k, names_prefix = "precast_k"
    )

temp <- colnames(spe)
colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(precast_results, by = "key") |>
    DataFrame()
colnames(spe) <- temp

################################################################################
#   Check that k = 2 strongly correlates with WM vs GM boundary
################################################################################

#   For each spot, average expression Z-scores across all white matter genes.
#   The idea is that the combined expression of several white matter genes
#   should form a more precise marker of white matter than one gene
a <- assays(spe)$logcounts[
    match(subregion_genes[["white_matter"]], rowData(spe)$gene_name),
]
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

################################################################################
#   Simply plot cluster assignments for k=2 through k=28
################################################################################

for (k in 2:28) {
    plot_list <- list()
    for (donor in unique(spe$sample_id)) {
        plot_list[[donor]] <- spot_plot(
            spe,
            sample_id = donor, var_name = paste0("precast_k", k),
            is_discrete = TRUE
        )
    }

    pdf(file.path(plot_dir, sprintf("k%s_all_samples.pdf", k)))
    print(plot_list)
    dev.off()
}

################################################################################
#   Evaluate if PRECAST cluster assignments were shared for overlapping spots
################################################################################

#   Given a character(1) 'this_key' (a unique barcode in 'col_data$key'), and
#   a tibble of colData 'col_data', return a character(1) giving the key of the
#   spot on a different capture area overlapping the spot given by 'this_key',
#   or return "" if no such spot exists
get_overlapping_key <- function(this_key, col_data) {
    i <- match(this_key, col_data$key)
    x <- col_data |>
        filter(
            sample_id_original == col_data[[i, "overlap_slide"]],
            array_row_transformed == col_data[[i, "array_row_transformed"]],
            array_col_transformed == col_data[[i, "array_col_transformed"]]
        ) |>
        #   Generally, we only had at most 2 capture areas overlapping, but due
        #   to rounding complexities there occasionally can be more than one
        #   spot for the same capture area and array coordinates (hence taking
        #   the first row)
        slice_head(n = 1) |>
        pull(key)

    if (length(x) == 0) {
        return("")
    } else {
        return(x)
    }
}

col_data_full <- colData(spe) |>
    as_tibble()

#   Here, just take the non-excluded spots that overlap another capture area
col_data_small <- col_data_full |>
    filter(!is.na(overlap_slide), overlap_slide != "None", !exclude_overlapping) |>
    select(matches("key|^sample_id$|^precast_k")) |>
    #   Effectively removes spots where PRECAST has not assigned a cluster
    #   identity
    na.omit()

#   Add the key of each overlapping spot, and remove spots where a key doesn't
#   exist
col_data_small <- col_data_small |>
    mutate(
        overlap_key = map_chr(
            col_data_small$key, get_overlapping_key,
            col_data = col_data_full
        )
    ) |>
    filter(overlap_key != "") |>
    #   Append '_original' to the PRECAST clustering results, to signify these
    #   are the results for the non-excluded spots
    rename_with(~ ifelse(grepl("^precast_k", .x), paste0(.x, "_original"), .x))

results_df <- col_data_small |>
    #   Add PRECAST results from the corresponding excluded overlapping spots,
    #   adding '_overlap' to the colnames to differentiate them
    cbind(
        col_data_full[match(col_data_small$overlap_key, col_data_full$key), ] |>
            select(matches("^precast_k")) |>
            rename_with(~ paste0(.x, "_overlap"))
    ) |>
    as_tibble() |>
    #   The colnames now include info about both the capture area (original vs.
    #   overlap) and value of k. Pivot longer to break into 3 columns: capture
    #   area, k, and the cluster identity (assignment)
    pivot_longer(
        cols = matches("^precast_k[0-9]+_"),
        names_to = c("k", "capture_area"),
        names_pattern = "^precast_k([0-9]+)_(overlap|original)$",
        values_to = "cluster_assignment"
    ) |>
    #   Pivot wider so we can compare the cluster identities for the original
    #   and overlapping spots side by side
    pivot_wider(
        names_from = "capture_area", values_from = "cluster_assignment"
    ) |>
    #   Fix a data type
    mutate(k = factor(as.integer(k))) |>
    #   Rarely, an overlapping spot but not the original may have been dropped
    #   as input to PRECAST; simply drop rows with this situation
    na.omit()

#   Plot the proportion of (unique) spots that match their overlapping spot.
#   Each boxplot contains all donors, and we split by value of k
p <- results_df |>
    group_by(sample_id, k) |>
    summarize(match_rate = mean(original == overlap)) |>
    ungroup() |>
    ggplot() +
    geom_boxplot(aes(x = k, y = match_rate), outlier.shape = NA) +
    labs(x = "PRECAST k Value", y = "Proportion of Matches", color = "Donor") +
    geom_jitter(aes(x = k, y = match_rate, color = sample_id)) +
    theme_bw(base_size = 24)

pdf(
    file.path(plot_dir, "cluster_assignment_at_overlaps.pdf"),
    width = 15, height = 7
)
print(p)
dev.off()

#   After some interactive exploration, the proportion of overlapping spots
#   belonging to white matter seems to be the best predictor I found of
#   cluster-assignment mismatches at overlap

#   First, prove the relationship exists at the spot level
message("Among overlapping spots, correlation between matching cluster assignment at")
message("k = 2 and at least one spot in the pair belonging to white matter:")
results_df |>
    filter(k == 2) |>
    summarize(
        wm_vs_match = cor(original == overlap, (original == 2) | (overlap == 2))
    ) |>
    pull(wm_vs_match) |>
    round(2) |>
    message()

#   Next, visually show how the effect manifests at the donor level
p = results_df |>
    filter(k == 2) |>
    group_by(sample_id) |>
    summarize(
        match_rate = mean(original == overlap),
        prop_white = mean((original == 2) | (overlap == 2))
    ) |>
    ungroup() |>
    ggplot() +
        geom_point(aes(x = prop_white, y = match_rate, color = sample_id)) +
        labs(
            x = 'Prop. Pairs w/ at Least 1 WM Spot',
            y = 'Prop. Pairs w/ Matching k = 2 Assignment',
            color = 'Donor'
        ) +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, "WM_vs_match_rate.pdf"))
print(p)
dev.off()

################################################################################
#   Compare PRECAST clusters with marker gene expression of suspected biological
#   regions (e.g. WM/GM, D1 islands, shell, core)
################################################################################

#   Specifically, k = 4 seems to consistently capture WM vs. GM, has a cluster
#   spread in small pockets across the tissue (D1 islands?), and a cluster
#   marking a subregion of GM (shell or core?)
#
#   Evaluate how well marker genes (albeit derived from mouse/ other species)
#   of these biological regions line up with k = 4 clusters from PRECAST
b <- Sys.time()
assay_subset <- Matrix(
    assays(spe)$logcounts[match(unlist(subregion_genes), rowData(spe)$gene_name), ],
    sparse = TRUE,
    nrow = length(unlist(subregion_genes)),
    dimnames = list(unlist(subregion_genes), colnames(spe))
)
a <- assay_subset[subregion_genes[[names(subregion_genes)[[1]]]], ]
gene_z <- colMeans((a - rowMeans(a)) / rowSds(a), na.rm = TRUE)
Sys.time() - b

#   Bring in the subset of the logcounts that we'll use, as this seems to be
#   significantly faster than repeatedly performing HDF5-backed, chunked
#   operations
assay_subset <- Matrix(
    assays(spe)$logcounts[match(unlist(subregion_genes), rowData(spe)$gene_name), ],
    sparse = TRUE,
    nrow = length(unlist(subregion_genes)),
    dimnames = list(unlist(subregion_genes), colnames(spe))
)

cor_mat <- matrix(0, nrow = length(subregion_genes), ncol = 4)
rownames(cor_mat) <- names(subregion_genes)
for (subregion in names(subregion_genes)) {
    #   For each spot, average expression Z-scores across the subregion.
    #   The idea is that the combined expression of several genes
    #   should form a more precise marker of the subregion than one gene
    a <- assay_subset[subregion_genes[[subregion]], ]
    gene_z <- colMeans((a - rowMeans(a)) / rowSds(a), na.rm = TRUE)

    #   Now compute correlation between cluster identity (in a boolean sense--
    #   belonging or not belonging in the cluster) and the averaged Z-score
    for (cluster_identity in 1:4) {
        cor_mat[match(subregion, names(subregion_genes)), cluster_identity] <-
            cor(
                gene_z, spe$precast_k4 == cluster_identity,
                use = "complete.obs"
            )
    }
}

message("Correlation of subregion markers with clusters in k=4:")
print(cor_mat)

session_info()
