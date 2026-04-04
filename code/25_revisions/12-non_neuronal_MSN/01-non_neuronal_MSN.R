rm(list = ls())
library(SpatialExperiment)
library(HDF5Array)
library(lobstr)
library(here)
library(sessioninfo)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(purrr)
library(MERINGUE)
library(spatialNAcUtils)
library(cowplot)
library(spatialLIBD)

spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)

plotDir <- here("plots", "25_revisions", "12-non_neuronal_MSN")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

spe <- spe[, !is.na(spe[["spatial_domains"]])]
spe <- spe[, spe$spatial_domains %in% c("MSN_1", "MSN_2", "MSN_3")]

plot_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank()
)

legend_theme <- guides(
  fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(2.5, "cm"),
    barwidth = unit(0.5, "cm")
  ),
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(2.5, "cm"),
    barwidth = unit(0.5, "cm")
  )
)

p_Oligo <- vis_gene(
  spe,
  sampleid = "Br6522",
  geneid = "Oligo",
  is_stitched = TRUE,
  cont_colors = viridisLite::rocket(10, direction = -1)
) +
  ggtitle("Oligo") +
  plot_theme +
  legend_theme

p_Astro_A <- vis_gene(
  spe,
  sampleid = "Br6522",
  geneid = "Astrocyte_A",
  is_stitched = TRUE,
  cont_colors = viridisLite::rocket(10, direction = -1)
) +
  ggtitle("Astrocyte A") +
  plot_theme +
  legend_theme

p_Astro_B <- vis_gene(
  spe,
  sampleid = "Br6522",
  geneid = "Astrocyte_B",
  is_stitched = TRUE,
  cont_colors = viridisLite::rocket(10, direction = -1)
) +
  ggtitle("Astrocyte B") +
  plot_theme +
  legend_theme

p_spot <- cowplot::plot_grid(
  p_Oligo, p_Astro_A, p_Astro_B,
  nrow = 1,
  align = "h"
)

pdf(file.path(plotDir, "RCTD_non_neuronal.pdf"), width = 12, height = 4)
print(p_spot)
dev.off()

load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)
celltype_colors <- cluster_cols[!names(cluster_cols) == "Neuro_Ambig"]

## -----------------------------------------
## Composition analysis (Astrocyte_A + Oligo)
## -----------------------------------------

df <- as.data.frame(colData(spe))

df <- df %>%
  filter(spatial_domains %in% c("MSN_1","MSN_2","MSN_3"))

msn_types <- c(
  "DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
  "DRD2_MSN_A","DRD2_MSN_B"
)

msn_colors <- celltype_colors[msn_types]

## -----------------------------
## Helper function
## -----------------------------

run_composition_plot <- function(df, celltype_name, threshold = 0.1, prefix) {

  df_sub <- df %>%
    filter(.data[[celltype_name]] > threshold)

  plot_df <- df_sub %>%
    select(spatial_domains, all_of(msn_types)) %>%
    pivot_longer(cols = all_of(msn_types),
                 names_to = "MSN",
                 values_to = "weight")

  plot_df$MSN <- factor(plot_df$MSN, levels = msn_types)

  ## pairwise stats within each spatial domain
  stat_df <- plot_df %>%
    group_by(spatial_domains) %>%
    do({
      d <- .
      combs <- combn(msn_types, 2, simplify = FALSE)

      res <- map_dfr(combs, function(pair) {
        x <- d$weight[d$MSN == pair[1]]
        y <- d$weight[d$MSN == pair[2]]

        if (length(x) < 5 | length(y) < 5) return(NULL)

        wt <- wilcox.test(x, y)

        data.frame(
          group1 = pair[1],
          group2 = pair[2],
          p = wt$p.value
        )
      })

      res$p_adj <- p.adjust(res$p, method = "BH")
      res
    }) %>%
    ungroup()

  ## plot
  p <- ggplot(plot_df, aes(x = MSN, y = weight, fill = MSN)) +
    geom_boxplot(outlier.size = 0, width = 0.7) +
    facet_wrap(~ spatial_domains, nrow = 1) +
    scale_fill_manual(values = msn_colors, drop = FALSE) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = paste0(prefix, " spots (>", threshold, ")"),
      x = NULL,
      y = "RCTD weight"
    )

  ggsave(
    file.path(plotDir, paste0(prefix, "_MSN_violin.pdf")),
    p,
    width = 8,
    height = 4
  )

  write.csv(
    stat_df,
    file.path(plotDir, paste0(prefix, "_pairwise_stats.csv")),
    row.names = FALSE
  )

  return(list(plot = p, stats = stat_df))
}

## -----------------------------
## Run analyses
## -----------------------------

astro_res <- run_composition_plot(
  df,
  celltype_name = "Astrocyte_A",
  threshold = 0.1,
  prefix = "AstrocyteA"
)

oligo_res <- run_composition_plot(
  df,
  celltype_name = "Oligo",
  threshold = 0.1,
  prefix = "Oligo"
)