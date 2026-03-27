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

# Read in SPE object
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)

plotDir <- here("plots", "25_revisions", "10-rebuttal_figure_RCTD")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

# Get the plot data
plot_df <- colData(spe)
plot_df <- plot_df[ ,c("spatial_domains", "Excitatory", "DRD1_MSN_A", "DRD2_MSN_A")]
plot_df <- plot_df[plot_df$spatial_domains %in% c("MSN_1", "MSN_2", "MSN_3"), ]

load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)

# Ensure this is a plain data.frame for tidyverse reshaping
plot_df <- as.data.frame(plot_df)

# Order MSN domains on the x-axis
plot_df$spatial_domains <- factor(
  plot_df$spatial_domains,
  levels = c("MSN_1", "MSN_2", "MSN_3")
)

# Reshape to long format
plot_long <- plot_df %>%
  pivot_longer(
    cols = c("Excitatory", "DRD1_MSN_A", "DRD2_MSN_A"),
    names_to = "cell_type",
    values_to = "weight"
  )

# Keep cell type order consistent
plot_long$cell_type <- factor(
  plot_long$cell_type,
  levels = c("Excitatory", "DRD1_MSN_A", "DRD2_MSN_A")
)

# Simple grouped boxplot
p <- ggplot(plot_long, aes(x = spatial_domains, y = weight, fill = cell_type)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = cluster_cols[c("Excitatory", "DRD1_MSN_A", "DRD2_MSN_A")]
  ) +
  labs(
    x = "MSN spatial domain",
    y = "RCTD weight",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

pdf(file.path(plotDir, "RCTD_weights.pdf"), width = 8, height = 2)
print(p)
dev.off()