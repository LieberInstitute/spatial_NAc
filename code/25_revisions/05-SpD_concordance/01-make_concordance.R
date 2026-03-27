library(here)
library(HDF5Array)
library(tidyverse)
library(SpatialExperiment)
library(spatialNAcUtils)
library(spatialLIBD)
library(ggpubr)
library(cowplot)
library(purrr)

## -----------------------------
## Paths
## -----------------------------
spe_dir <- here(
    "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)

precast_dir <- here(
    "processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast"
)

plot_dir <- here(
    "plots", "25_revisions", "05-SpD_concordance"
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## Load SPE
## -----------------------------
spe <- loadHDF5SummarizedExperiment(spe_dir)

## -----------------------------
## Read PRECAST k = 10 results
## -----------------------------
k <- 10
n_random_starts <- 5

precast_list <- lapply(seq_len(n_random_starts), function(rs) {
    fp <- file.path(precast_dir, paste0("random_start_", rs), sprintf("PRECAST_k%s.csv", k))

    read.csv(fp) |>
        as_tibble() |>
        select(key, cluster) |>
        mutate(cluster = as.character(cluster)) |>
        rename(!!paste0("cluster_rs", rs) := cluster)
})

precast_df <- purrr::reduce(precast_list, left_join, by = "key")

## -----------------------------
## Check keys are consistent
## -----------------------------
stopifnot(all(map_lgl(precast_list[-1], ~ identical(precast_list[[1]]$key, .x$key))))

## -----------------------------
## Domain mapping by random start
## Conservative mapping:
## - unassigned / ambiguous clusters remain NA
## - random start 3 cluster 3 is left unassigned
## -----------------------------
domain_map <- list(
    rs1 = c(
        "3" = "MSN_2/MSN_1",
        "4" = "MSN_3",
        "8" = "D1_islands"
    ),
    rs2 = c(
        "5" = "MSN_1",
        "1" = "MSN_2",
        "3" = "MSN_3",
        "9" = "D1_islands"
    ),
    rs3 = c(
        "1" = "MSN_1",
        "4" = "MSN_2",
        "5" = "MSN_3",
        "3" = "D1_islands"
    ),
    rs4 = c(
        "1" = "MSN_1",
        "4" = "MSN_2",
        "2" = "MSN_3",
        "3" = "D1_islands"
    ),
    rs5 = c(
        "4" = "MSN_1",
        "3" = "MSN_2/MSN_3",
        "6" = "MSN_2",
        "9" = "MSN_3",
        "1" = "D1_islands"
    )
)

annotate_domain <- function(cluster_vec, mapping_vec) {
    out <- unname(mapping_vec[cluster_vec])
    out[is.na(out)] <- NA_character_
    out
}

## Annotate each random start
annotated_df <- precast_df |>
    mutate(
        domain_rs1 = annotate_domain(cluster_rs1, domain_map$rs1),
        domain_rs2 = annotate_domain(cluster_rs2, domain_map$rs2),
        domain_rs3 = annotate_domain(cluster_rs3, domain_map$rs3),
        domain_rs4 = annotate_domain(cluster_rs4, domain_map$rs4),
        domain_rs5 = annotate_domain(cluster_rs5, domain_map$rs5)
    )

## -----------------------------
## Compute per-spot proportions across 5 random starts
## -----------------------------
domain_cols <- paste0("domain_rs", 1:5)

calc_prop <- function(df, target_domain) {
    m <- as.matrix(df[, domain_cols])
    rowMeans(m == target_domain, na.rm = TRUE)
}

annotated_df$prop_MSN_1      <- calc_prop(annotated_df, "MSN_1")
annotated_df$prop_MSN_2      <- calc_prop(annotated_df, "MSN_2")
annotated_df$prop_MSN_3      <- calc_prop(annotated_df, "MSN_3")
annotated_df$prop_D1_islands <- calc_prop(annotated_df, "D1_islands")

## -----------------------------
## Add annotations + proportions to spe
## -----------------------------
temp <- colnames(spe)

colData(spe) <- colData(spe) |>
    as_tibble() |>
    left_join(annotated_df, by = "key") |>
    DataFrame()

colnames(spe) <- temp

## Also make a tibble version for plotting
cd <- as_tibble(colData(spe))

## -----------------------------
## 1. Density plots
## -----------------------------
density_df <- cd |>
    select(
        key,
        prop_MSN_1,
        prop_MSN_2,
        prop_MSN_3,
        prop_D1_islands
    ) |>
    pivot_longer(
        cols = starts_with("prop_"),
        names_to = "domain",
        values_to = "proportion"
    ) |>
    mutate(
        domain = recode(
            domain,
            prop_MSN_1 = "MSN_1",
            prop_MSN_2 = "MSN_2",
            prop_MSN_3 = "MSN_3",
            prop_D1_islands = "D1 islands"
        ),
        domain = factor(domain, levels = c("D1 islands", "MSN_1", "MSN_2", "MSN_3"))
    )

p_density <- ggplot(
    density_df %>% filter(proportion > 0),
    aes(x = proportion)
) +
    geom_histogram(
        binwidth = 0.1,   # thinner bars (0.1 instead of 0.2)
        boundary = 0,
        color = "black",
        fill = "grey40"
    ) +
    facet_wrap(~ domain, nrow = 1, scales = "free_y") +
    theme_pubr() +
    xlab("Proportion of random starts assigned to domain") +
    ylab("Count")

ggsave(
    filename = file.path(plot_dir, "density_domain_assignment_proportions_k10.pdf"),
    plot = p_density,
    width = 12,
    height = 3
)

# Add spot plots
p_D1_islands_spot <- vis_gene(spe, sampleid = "Br6522", geneid = "prop_D1_islands", is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
p_MSN_1_spot <- vis_gene(spe, sampleid = "Br6522", geneid = "prop_MSN_1", is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
p_MSN_2_spot <- vis_gene(spe, sampleid = "Br6522", geneid = "prop_MSN_2", is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))
p_MSN_3_spot <- vis_gene(spe, sampleid = "Br6522", geneid = "prop_MSN_3", is_stitched = TRUE, cont_colors = viridisLite::rocket(10, direction = -1))

pdf(file.path(plot_dir, "D1_spot_plot.pdf"), width = 4, height = 3)
print(p_D1_islands_spot)
dev.off()

pdf(file.path(plot_dir, "MSN_1_spot_plot.pdf"), width = 5, height = 4)
print(p_MSN_1_spot)
dev.off()

pdf(file.path(plot_dir, "MSN_2_spot_plot.pdf"), width = 5, height = 4)
print(p_MSN_2_spot)
dev.off()

pdf(file.path(plot_dir, "MSN_3_spot_plot.pdf"), width = 5, height = 4)
print(p_MSN_3_spot)
dev.off()