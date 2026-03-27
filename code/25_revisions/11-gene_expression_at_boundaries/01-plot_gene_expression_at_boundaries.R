rm(list = ls())
library("SpatialExperiment")
library("HDF5Array")
library("lobstr")
library("here")
library("sessioninfo")
library("Matrix")
library("dplyr")
library("ggplot2")
library("tidyr")
library("scales")
library("purrr")
library(MERINGUE)
library(spatialNAcUtils)

# =========================================================
# Load data
# =========================================================
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)
plotDir <- here("plots", "25_revisions", "11-gene_expression_at_boundaries")

spe <- spe[ ,!is.na(spe[["spatial_domains"]])]

# Keep only MSN domains of interest
spe <- spe[, spe$spatial_domains %in% c("MSN_1", "MSN_2", "MSN_3")]
spe <- spe[ ,!spe$exclude_overlapping]
rownames(spatialCoords(spe)) <- colnames(spe)

# =========================================================
# Run per donor and count neighbor domains
# =========================================================
sample_ids <- unique(as.character(spe$sample_id))

neighbor_count_list <- lapply(sample_ids, function(sample_id) {
    message("Processing sample: ", sample_id)

    spe_sub <- spe[, spe$sample_id == sample_id]

    W <- getSpatialNeighbors(spatialCoords(spe_sub), filterDist = 500)

    if (is.null(rownames(W))) {
        rownames(W) <- colnames(spe_sub)
        colnames(W) <- colnames(spe_sub)
    }

    spot_ids <- colnames(spe_sub)
    spot_domains <- spe_sub$spatial_domains
    names(spot_domains) <- spot_ids

    # count neighbors excluding self if present
    diag_present <- any(diag(W) != 0)

    sample_df <- lapply(seq_along(spot_ids), function(i) {
        spot_id <- spot_ids[i]

        nbr_idx <- which(W[spot_id, ] != 0)
        nbr_ids <- colnames(W)[nbr_idx]

        if (diag_present) {
            nbr_ids <- setdiff(nbr_ids, spot_id)
        }

        nbr_domains <- spot_domains[nbr_ids]

        data.frame(
            spot_id = spot_id,
            sample_id = sample_id,
            spot_domain = spot_domains[spot_id],
            n_neighbors = length(nbr_ids),
            n_nbr_MSN_1 = sum(nbr_domains == "MSN_1"),
            n_nbr_MSN_2 = sum(nbr_domains == "MSN_2"),
            n_nbr_MSN_3 = sum(nbr_domains == "MSN_3"),
            stringsAsFactors = FALSE
        )
    }) |> bind_rows()

    stopifnot(all(
        sample_df$n_neighbors ==
            sample_df$n_nbr_MSN_1 +
            sample_df$n_nbr_MSN_2 +
            sample_df$n_nbr_MSN_3
    ))

    sample_df
})

neighbor_count_df <- bind_rows(neighbor_count_list)

# quick overall check
head(neighbor_count_df)
summary(neighbor_count_df$n_neighbors)

# =========================================================
# Add back to spe
# =========================================================
match_idx <- match(colnames(spe), neighbor_count_df$spot_id)
stopifnot(!any(is.na(match_idx)))

spe$n_neighbors <- neighbor_count_df$n_neighbors[match_idx]
spe$n_nbr_MSN_1 <- neighbor_count_df$n_nbr_MSN_1[match_idx]
spe$n_nbr_MSN_2 <- neighbor_count_df$n_nbr_MSN_2[match_idx]
spe$n_nbr_MSN_3 <- neighbor_count_df$n_nbr_MSN_3[match_idx]

# =========================================================
# Define boundary labels
# =========================================================
spe$boundary_type <- NA_character_

# MSN_1 spots
idx <- spe$spatial_domains == "MSN_1"
spe$boundary_type[idx & spe$n_nbr_MSN_2 == 0 & spe$n_nbr_MSN_3 == 0] <- "MSN_1_interior"
spe$boundary_type[idx & spe$n_nbr_MSN_2 > 0 & spe$n_nbr_MSN_3 == 0] <- "MSN_1_bordering_MSN_2"
spe$boundary_type[idx & spe$n_nbr_MSN_2 == 0 & spe$n_nbr_MSN_3 > 0] <- "MSN_1_bordering_MSN_3"
spe$boundary_type[idx & spe$n_nbr_MSN_2 > 0 & spe$n_nbr_MSN_3 > 0] <- "MSN_1_bordering_MSN_2_MSN_3"

# MSN_2 spots
idx <- spe$spatial_domains == "MSN_2"
spe$boundary_type[idx & spe$n_nbr_MSN_1 == 0 & spe$n_nbr_MSN_3 == 0] <- "MSN_2_interior"
spe$boundary_type[idx & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_3 == 0] <- "MSN_2_bordering_MSN_1"
spe$boundary_type[idx & spe$n_nbr_MSN_1 == 0 & spe$n_nbr_MSN_3 > 0] <- "MSN_2_bordering_MSN_3"
spe$boundary_type[idx & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_3 > 0] <- "MSN_2_bordering_MSN_1_MSN_3"

# MSN_3 spots
idx <- spe$spatial_domains == "MSN_3"
spe$boundary_type[idx & spe$n_nbr_MSN_1 == 0 & spe$n_nbr_MSN_2 == 0] <- "MSN_3_interior"
spe$boundary_type[idx & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_2 == 0] <- "MSN_3_bordering_MSN_1"
spe$boundary_type[idx & spe$n_nbr_MSN_1 == 0 & spe$n_nbr_MSN_2 > 0] <- "MSN_3_bordering_MSN_2"
spe$boundary_type[idx & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_2 > 0] <- "MSN_3_bordering_MSN_1_MSN_2"

table(spe$boundary_type, useNA = "ifany")

# Also add pairwise boundaries

# =========================================================
# Pairwise boundary labels
# =========================================================
spe$pairwise_boundary <- NA_character_

# MSN_1 <-> MSN_2
spe$pairwise_boundary[
    (spe$spatial_domains == "MSN_1" & spe$n_nbr_MSN_2 > 0 & spe$n_nbr_MSN_3 == 0) |
    (spe$spatial_domains == "MSN_2" & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_3 == 0)
] <- "MSN_1_MSN_2"

# MSN_1 <-> MSN_3
spe$pairwise_boundary[
    (spe$spatial_domains == "MSN_1" & spe$n_nbr_MSN_3 > 0 & spe$n_nbr_MSN_2 == 0) |
    (spe$spatial_domains == "MSN_3" & spe$n_nbr_MSN_1 > 0 & spe$n_nbr_MSN_2 == 0)
] <- "MSN_1_MSN_3"

# MSN_2 <-> MSN_3
spe$pairwise_boundary[
    (spe$spatial_domains == "MSN_2" & spe$n_nbr_MSN_3 > 0 & spe$n_nbr_MSN_1 == 0) |
    (spe$spatial_domains == "MSN_3" & spe$n_nbr_MSN_2 > 0 & spe$n_nbr_MSN_1 == 0)
] <- "MSN_2_MSN_3"

# spots bordering both other domains
spe$pairwise_boundary[
    is.na(spe$pairwise_boundary) & grepl("bordering_.*_.*", spe$boundary_type)
] <- "multi_boundary"

# interior spots
spe$pairwise_boundary[grepl("_interior$", spe$boundary_type)] <- "interior"
spe$pairwise_boundary_w_domains <- spe$pairwise_boundary
spe$pairwise_boundary_w_domains[spe$pairwise_boundary == "interior" & spe$spatial_domains == "MSN_1"] <- "MSN_1"
spe$pairwise_boundary_w_domains[spe$pairwise_boundary == "interior" & spe$spatial_domains == "MSN_2"] <- "MSN_2"
spe$pairwise_boundary_w_domains[spe$pairwise_boundary == "interior" & spe$spatial_domains == "MSN_3"] <- "MSN_3"

table(spe$pairwise_boundary_w_domains, useNA = "ifany")
spe <- spe[ ,spe$n_neighbors >= 3]

# Check spot plots
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

sample_order <- unique(as.character(spe$sample_id))

# set factor order
plot_levels <- c(
    "MSN_1",
    "MSN_2",
    "MSN_3",
    "MSN_1_MSN_2",
    "MSN_1_MSN_3",
    "MSN_2_MSN_3",
    "multi_boundary"
)
plot_levels <- plot_levels[plot_levels %in% unique(spe$pairwise_boundary_w_domains)]

spe$pairwise_boundary_w_domains <- factor(
    spe$pairwise_boundary_w_domains,
    levels = plot_levels
)

# colors
plot_cols <- c(
    "MSN_1" = "#1b9e77",
    "MSN_2" = "#d95f02",
    "MSN_3" = "#7570b3",
    "MSN_1_MSN_2" = "#66a61e",
    "MSN_1_MSN_3" = "#e7298a",
    "MSN_2_MSN_3" = "#e6ab02",
    "multi_boundary" = "#a6761d"
)
plot_cols <- plot_cols[names(plot_cols) %in% plot_levels]

for (donor in sample_order) {

    message("Plotting donor: ", donor)

    p <- spot_plot(
        spe,
        sample_id = donor,
        var_name = "pairwise_boundary_w_domains",
        is_discrete = TRUE,
        spatial = TRUE,
        colors = plot_cols
    ) +
        ggtitle(donor) +
        guides(fill = guide_legend(override.aes = list(size = 5))) +
        theme(
            plot.title = element_text(face = "bold"),
            legend.position = "right"
        )

    # set size conditionally
    if (donor == "Br2720") {
        w <- 8
        h <- 8
    } else {
        w <- 9
        h <- 9
    }

    pdf(
        file.path(plotDir, paste0("pairwise_boundary_", donor, ".pdf")),
        width = w, height = h
    )
    print(p)
    dev.off()
}

# overall counts
count_df <- data.frame(
    pairwise_boundary_w_domains = spe$pairwise_boundary_w_domains
) |>
    dplyr::count(pairwise_boundary_w_domains, name = "n_spots") |>
    dplyr::filter(!is.na(pairwise_boundary_w_domains))

p_bar <- ggplot(count_df, aes(x = pairwise_boundary_w_domains, y = n_spots, fill = pairwise_boundary_w_domains)) +
    geom_col() +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = plot_cols[names(plot_cols) %in% count_df$pairwise_boundary_w_domains]) +
    labs(
        x = NULL,
        y = "Number of spots",
        title = "Spots per pairwise boundary category"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(face = "bold")
    )

ggsave(
    filename = file.path(plotDir, "barplot_pairwise_boundary_w_domains_counts.pdf"),
    plot = p_bar,
    width = 7,
    height = 4.5
)


# =========================================================
# Multi-panel violin plots for boundary categories
# =========================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

genes <- c("DRD1", "DRD2", "CALB1", "CARTPT", "TAC1", "PDYN", "ADORA2A", "PDE10A")

plot_levels <- c(
    "MSN_1",
    "MSN_1_MSN_2",
    "MSN_2",
    "MSN_2_MSN_3",
    "MSN_3",
    "MSN_1_MSN_3",
    "multi_boundary"
)
plot_levels <- plot_levels[plot_levels %in% unique(as.character(spe$pairwise_boundary_w_domains))]

spe$pairwise_boundary_w_domains <- factor(
    as.character(spe$pairwise_boundary_w_domains),
    levels = plot_levels
)

plot_cols <- c(
    "MSN_1" = "#1b9e77",
    "MSN_1_MSN_2" = "#66a61e",
    "MSN_2" = "#d95f02",
    "MSN_2_MSN_3" = "#e6ab02",
    "MSN_3" = "#7570b3",
    "MSN_1_MSN_3" = "#e7298a",
    "multi_boundary" = "#a6761d"
)
plot_cols <- plot_cols[names(plot_cols) %in% plot_levels]

# make long dataframe
expr_mat <- logcounts(spe)[genes, , drop = FALSE]

plot_df <- as.data.frame(t(as.matrix(expr_mat)))
colnames(plot_df) <- genes
plot_df$pairwise_boundary_w_domains <- spe$pairwise_boundary_w_domains

plot_df <- plot_df |>
    tidyr::pivot_longer(
        cols = all_of(genes),
        names_to = "gene",
        values_to = "logcounts"
    ) |>
    dplyr::filter(!is.na(pairwise_boundary_w_domains))

# optional cleaner labels
plot_df$pairwise_boundary_w_domains <- factor(
    plot_df$pairwise_boundary_w_domains,
    levels = plot_levels
)

make_violin_plot <- function(gene_name, panel_label = NULL) {
    df_sub <- plot_df |> dplyr::filter(gene == gene_name)

    p <- ggplot(df_sub, aes(x = pairwise_boundary_w_domains, y = logcounts, fill = pairwise_boundary_w_domains)) +
        geom_violin(scale = "width", trim = TRUE, color = "black", size = 0.3) +
        scale_fill_manual(values = plot_cols, drop = FALSE) +
        theme_classic(base_size = 14) +
        labs(
            title = gene_name,
            x = NULL,
            y = "logcounts"
        ) +
        theme(
            plot.title = element_text(face = "italic", hjust = 0.5, size = 18),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 16),
            legend.position = "none"
        )

    if (!is.null(panel_label)) {
        p <- p + labs(subtitle = panel_label) +
            theme(
                plot.subtitle = element_text(face = "bold", size = 20, hjust = -0.05)
            )
    }

    p
}

panel_labels <- LETTERS[seq_along(genes)+2]

plot_list <- purrr::map2(
    genes, panel_labels,
    ~make_violin_plot(.x, .y)
)

p_combined <- wrap_plots(plot_list, ncol = 2)

ggsave(
    filename = file.path(plotDir, "violin_plots_pairwise_boundary_w_domains.pdf"),
    plot = p_combined,
    width = 14,
    height = 18
)