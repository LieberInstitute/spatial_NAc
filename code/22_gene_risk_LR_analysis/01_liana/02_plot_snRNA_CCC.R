library(liana)
library(tidyverse)
library(unglue)
library(SingleCellExperiment)
library(here)
library(ComplexHeatmap)
library(circlize)
library(glue)

plot_liana_heatmap <- function(mat, font_size = 12, grid_text = FALSE,
    name = "Frequency",
    pallette = c("white", "violetred2"),
    row_title = "Sender (Cell types)",
    column_title = "Receiver (Cell types)",
    cell_col = NULL, # Boyi Added parameter
    ...) {
    if (grid_text) {
        grid_text <- function(j, i, x, y, width, height, fill) {
            grid_text <- grid.text(sprintf("%d", mat[i, j]),
                x, y,
                gp = gpar(fontsize = font_size * 0.83)
            )
        }
    } else {
        grid_text <- NULL
    }

    # browser()
    cell_anno <- unique(rownames(mat))
    if (is.null(cell_col)) {
        fun_color <- grDevices::colorRampPalette(
            RColorBrewer::brewer.pal(n = 8, name = "Dark2")
        )
        cell_anno <- fun_color(length(cell_anno)) %>%
            setNames(cell_anno)
    } else {
        cell_anno <- cell_col
    }

    ha_opts <- list(
        show_legend = FALSE,
        show_annotation_name = FALSE,
        col = list(anno = cell_anno),
        simple_anno_size = grid::unit(0.25, "cm")
    )

    column_ha <- exec("HeatmapAnnotation",
        anno = names(cell_anno),
        !!!ha_opts
    )

    row_ha <- exec("rowAnnotation",
        anno = names(cell_anno),
        !!!ha_opts
    )

    column_bar <- ComplexHeatmap::HeatmapAnnotation(
        bar = liana:::.anno_barplot(colSums(mat),
            cell_anno,
            axis.font.size = font_size * 0.4
        ),
        annotation_name_gp = gpar(fontsize = font_size * 0.5),
        show_legend = FALSE, show_annotation_name = FALSE
    )

    row_bar <- ComplexHeatmap::rowAnnotation(
        bar2 = liana:::.anno_barplot(
            rowSums(mat),
            cell_anno, font_size * 0.4
        ),
        gp = gpar(fill = cell_anno, col = cell_anno),
        show_legend = FALSE, show_annotation_name = FALSE # ,
        # gap = unit(10,"points")
    )

    ComplexHeatmap::Heatmap(
        mat,
        col = colorRampPalette(pallette)(10),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_side = "left",
        top_annotation = column_bar,
        bottom_annotation = column_ha,
        right_annotation = row_bar,
        left_annotation = row_ha,
        row_title = row_title,
        row_names_gp = gpar(fontsize = font_size),
        row_title_gp = gpar(fontsize = font_size * 1.2),
        column_names_gp = gpar(fontsize = font_size),
        column_title = column_title,
        column_title_gp = gpar(fontsize = font_size * 1.2),
        column_title_side = "bottom",
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = font_size *
                0.9, fontface = "bold"),
            border = NA, labels_gp = gpar(fontsize = font_size * 0.9),
            grid_width = unit(2, "mm")
        ),
        name = name, cell_fun = grid_text,
        ...
    )
}

plot_liana_circ <- function(liana_res, source_groups = NULL, target_groups = NULL,
    cex = 2, transparency = 0.4,
    facing = "clockwise", adj = c(-0.5, 0.05),
    cell_col = NULL, ...) {
    freqs <- liana_res %>%
        liana:::.get_freq()
    celltypes <- union(colnames(freqs), rownames(freqs))
    if (is.null(cell_col)) {
        grid.col <- (grDevices::colorRampPalette(
            (RColorBrewer::brewer.pal(n = 8, name = "Dark2"))
        )
        )(length(celltypes)) %>%
            setNames(celltypes)
    } else {
        grid.col <- cell_col
    }


    circlize::circos.clear()
    circlize::circos.par(circle.margin = 1)
    circlize::chordDiagram(freqs,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        link.arr.type = "big.arrow",
        transparency = transparency,
        grid.col = grid.col,
        annotationTrack = c("grid"), self.link = 1,
        big.gap = 7.5, small.gap = 5, ...
    )
    circlize::circos.trackPlotRegion(
        track.index = 1,
        panel.fun = function(x, y) {
            xlim <- circlize::get.cell.meta.data("xlim")
            ylim <- circlize::get.cell.meta.data("ylim")
            sector.name <- circlize::get.cell.meta.data("sector.index")
            circlize::circos.text(mean(xlim), ylim[1], sector.name,
                facing = facing, niceFacing = TRUE, adj = adj, cex = cex
            )
        }, bg.border = NA
    )
    p <- grDevices::recordPlot()
    return(p)
}

# Read in the single cell data
dat_dir <- here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(dat_dir, "sce_CellType_noresiduals.Rds"))

# Read cell type colors
cell_type_colors <-  c(
    "Astrocyte_A" = "#D95F00",
    "Astrocyte_B" = "#FF0016",
    "DRD1_MSN_A" = "#ECA31C",
    "DRD1_MSN_B" = "#D079AA",
    "DRD1_MSN_C" = "#352391",
    "DRD1_MSN_D" = "#d156aa",
    "DRD2_MSN_A" = "#58B6ED",
    "DRD2_MSN_B" ="#F80091",
    "Endothelial" = "#D00DFF",
    "Ependymal" = "#0D73B4",
    "Excitatory" = "#5C6300",
    "Inh_A" = "#35FB00",
    "Inh_B" = "#7A0096",
    "Inh_C" = "#854222", 
    "Inh_D" = "#A7F281", 
    "Inh_E" = "#0DFBFA", 
    "Inh_F" = "black",
    "Microglia" = "#F2E642",
    "Oligo" = "#4F4753",
    "OPC" = "#0D9F72"
)

plot_path <- here("plots", "22_gene_risk_LR_analysis", "01_liana")
res_dir <- here("processed-data", "22_gene_risk_LR_analysis", "01_liana")
liana_res <- readRDS(file.path(res_dir, "liana_consensus.rds"))

liana_trunc <- liana_res %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01)

liana_trunc |>
    arrange(aggregate_rank) |>
    write.csv(
        file = here(
            "processed-data/22_gene_risk_LR_analysis/01_liana/",
            "Table_LIANA_LR.csv"
        ),
        row.names = FALSE
    )

nrow(liana_trunc)

liana_trunc |>
    transmute(
        pair = glue::glue("{ligand.complex}_{receptor.complex}")
    ) |>
    summarize(n_distinct(pair))

pdf(
    file = file.path(plot_path, "liana_heatmap_overall.pdf"),
    width = 6,
    height = 6
)
plot_liana_heatmap(
    liana_trunc |>
        mutate(
            source = factor(
                source,
                levels = c(
                    "Astrocyte_A", "Astrocyte_B", "DRD1_MSN_A",
                    "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D",
                    "DRD2_MSN_A", "DRD2_MSN_B", "Excitatory", 
                    "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E",
                    "Inh_F", "Microglia", "Oligo", "OPC", "Endothelial", 
                    "Ependymal"
                )
            ),
            target = factor(
                target,
                levels = c(
                   "Astrocyte_A", "Astrocyte_B", "DRD1_MSN_A",
                    "DRD1_MSN_B", "DRD1_MSN_C", "DRD1_MSN_D",
                    "DRD2_MSN_A", "DRD2_MSN_B", "Excitatory", 
                    "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", 
                    "Inh_F", "Microglia", "Oligo", "OPC", "Endothelial", 
                    "Ependymal"
                )
            )
        ) |>
        liana:::.get_freq(),
    cell_col = cell_type_colors # ,
    # row_dend_reorder = FALSE,
    # column_dend_reorder = FALSE
) |>
    print()
dev.off()

# Ligands of interest
# Define ligands of interest
ligands_of_interest <- c("SEMA3E", "TAC1", "PENK")

# Count interactions by ligand and source (sender) cell type
ligand_source_counts <- liana_trunc |>
  dplyr::filter(ligand.complex %in% ligands_of_interest) |>
  dplyr::count(ligand.complex, source) |>
  dplyr::filter(n > 0)  # Remove truly zero counts

# Plot with free x- and y-axis scaling per facet
p <- ggplot(ligand_source_counts, aes(x = source, y = n, fill = source)) +
  geom_col() +
  facet_wrap(~ ligand.complex, scales = "free") +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 1)
  ) +
  scale_fill_manual(values = cell_type_colors) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  labs(
    x = "Sender Cell Type",
    y = "Interaction Count",
    title = "Sender Cell Type Frequency per Ligand"
  )
# Save it
ggsave(
  filename = file.path(plot_path, "grouped_barplot_sender_per_ligand.pdf"),
  plot = p,
  width = 9,
  height = 3.5
)

# Receptors of interest
# Define receptors of interest
receptors_of_interest <- c("OPRM1", "PLXND1", "TACR3")

# Count interactions by ligand and source (sender) cell type
receptor_target_counts <- liana_trunc |>
  dplyr::filter(receptor.complex %in% receptors_of_interest) |>
  dplyr::count(receptor.complex, target) |>
  dplyr::filter(n > 0)  # Remove truly zero counts

# Plot with free x- and y-axis scaling per facet
p <- ggplot(receptor_target_counts, aes(x = target, y = n, fill = target)) +
  geom_col() +
  facet_wrap(~ receptor.complex, scales = "free") +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 1)
  ) +
  scale_fill_manual(values = cell_type_colors) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  ) +
  labs(
    x = "Receiver Cell Type",
    y = "Interaction Count",
    title = "Receiver Cell Type Frequency per Receptor"
  )
# Save it
ggsave(
  filename = file.path(plot_path, "grouped_barplot_target_per_receptor.pdf"),
  plot = p,
  width = 9,
  height = 3.5
)



# Filter for NPY interactions
npy_res <- liana_trunc %>%
  filter(ligand.complex == "NPY")
npy_res$LR_pair <- paste0(npy_res$ligand.complex, " -> ", npy_res$receptor.complex)

# Generate dotplot with top 10 interactions
p <- ggplot(npy_res, aes(
        x = target,
        y = LR_pair,
        size = natmi.edge_specificity,
        color = sca.LRscore
    )) +
        geom_point() +
        scale_color_viridis_c(name = "Expression\nMagnitude", option = "viridis") +
        scale_size(range = c(2, 6), name = "Interaction\nSpecificity") +
        labs(
            x = "Receiver (Target)",
            y = ""
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(drop = TRUE)

# Save the plot
ggsave(
  filename = here::here("plots", "22_gene_risk_LR_analysis", "01_liana", "dotplot_NPY_InhE.pdf"),
  plot = p,
  width = 7,
  height = 5
)

tac1_res <- liana_trunc %>%
  filter(ligand.complex == "TAC1") %>%
  filter(source == "DRD2_MSN_B")
tac1_res$LR_pair <- paste0(tac1_res$ligand.complex, " -> ", tac1_res$receptor.complex)

# Generate dotplot with top 10 interactions
p <- ggplot(tac1_res, aes(
        x = target,
        y = LR_pair,
        size = natmi.edge_specificity,
        color = sca.LRscore
    )) +
        geom_point() +
        scale_color_viridis_c(name = "Expression\nMagnitude", option = "viridis") +
        scale_size(range = c(2, 6), name = "Interaction\nSpecificity") +
        labs(
            x = "Receiver (Target)",
            y = ""
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(drop = TRUE)

# Save the plot
ggsave(
  filename = here::here("plots", "22_gene_risk_LR_analysis", "01_liana", "dotplot_TAC1_DRD2_MSN_B.pdf"),
  plot = p,
  width = 7,
  height = 5
)

sema3e_res <- liana_trunc %>%
  filter(ligand.complex == "SEMA3E") %>%
  filter(source == "DRD1_MSN_C")
sema3e_res$LR_pair <- paste0(sema3e_res$ligand.complex, " -> ", sema3e_res$receptor.complex)

# Generate dotplot with top 10 interactions
p <- ggplot(sema3e_res, aes(
        x = target,
        y = LR_pair,
        size = natmi.edge_specificity,
        color = sca.LRscore
    )) +
        geom_point() +
        scale_color_viridis_c(name = "Expression\nMagnitude", option = "viridis") +
        scale_size(range = c(2, 6), name = "Interaction\nSpecificity") +
        labs(
            x = "Receiver (Target)",
            y = ""
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(drop = TRUE)

# Save the plot
ggsave(
  filename = here::here("plots", "22_gene_risk_LR_analysis", "01_liana", "dotplot_SEMA3E_DRD1_MSN_C.pdf"),
  plot = p,
  width = 7,
  height = 5
)

penk_res <- liana_trunc %>%
  filter(ligand.complex == "PENK") %>%
  filter(source == "DRD2_MSN_A")
penk_res$LR_pair <- paste0(penk_res$ligand.complex, " -> ", penk_res$receptor.complex)

# Generate dotplot with top 10 interactions
p <- ggplot(penk_res, aes(
        x = target,
        y = LR_pair,
        size = natmi.edge_specificity,
        color = sca.LRscore
    )) +
        geom_point() +
        scale_color_viridis_c(name = "Expression\nMagnitude", option = "viridis") +
        scale_size(range = c(2, 6), name = "Interaction\nSpecificity") +
        labs(
            x = "Receiver (Target)",
            y = ""
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(drop = TRUE)

# Save the plot
ggsave(
  filename = here::here("plots", "22_gene_risk_LR_analysis", "01_liana", "dotplot_PENK_DRD2_MSN_A.pdf"),
  plot = p,
  width = 7,
  height = 5
)

lr_pairs <- list(
    c("SEMA3E", "PLXND1"),
    c("EFNA5", "EPHA5"),
    c("PENK", "OPRM1"), 
    c("PDYN", "OPRM1"),
    c("CRH", "CRHR1"), 
    c("TAC1", "TACR3"), 
    c("BDNF", "NTRK2"), 
    c("NPY", "MC4R"), 
    c("NPY", "DPP4"),
    c("NLGN3", "NRXN1"), 
    c("NLGN1", "NRXN1"),
    c("NLGN1", "NRXN3"), 
    c("CNTN3", "PTPRG"), 
    c("CNTN4", "PTPRG"), 
    c("CNTN6", "PTPRG")
)

for (pair in lr_pairs) {
    lig <- pair[1]
    rec <- pair[2]

    # Subset the liana_trunc data
    df <- liana_trunc %>%
        filter(
            ligand.complex == lig,
            receptor.complex == rec
        ) %>%
        mutate(
            interaction = paste(ligand.complex, " -> ", receptor.complex),
            source = factor(source, levels = unique(source)),
            target = factor(target, levels = unique(target))
        )

    # Skip empty subsets
    if (nrow(df) == 0) {
        message("Skipping: No interactions found for ", lig, " - ", rec)
        next
    }

    # Count unique sender and receiver cell types for sizing
    n_sources <- length(unique(df$source))
    n_targets <- length(unique(df$target))
    
    # Dynamically adjust plot size
    plot_width <- 3.5 + 0.5 * n_targets
    plot_height <- 3.5 + 0.5 * n_sources

    # Generate the plot
    p <- ggplot(df, aes(
        x = target,
        y = source,
        size = natmi.edge_specificity,
        color = sca.LRscore
    )) +
        geom_point() +
        scale_color_viridis_c(name = "Expression\nMagnitude", option = "viridis") +
        scale_size(range = c(2, 10), name = "Interaction\nSpecificity") +
        labs(
            x = "Receiver (Target)",
            y = "Sender (Source)",
            title = unique(df$interaction)
        ) +
        theme_bw(base_size = 14) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_y_discrete(drop = TRUE) +
        scale_x_discrete(drop = TRUE)

    # Save to PDF with dynamic size
    ggsave(
        filename = file.path(plot_path, paste0("dotplot_", lig, "_", rec, ".pdf")),
        plot = p,
        width = plot_width,
        height = plot_height
    )
}
