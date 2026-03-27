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

spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)

plot_sample_order <- c(
  "Br2743", "Br6432", "Br6423", "Br2720", "Br6471",
  "Br6522", "Br8492", "Br8325", "Br8667", "Br3942"
)

domain_levels <- c(
  "D1 islands", "Endothelial/Ependymal", "Excitatory",
  "Inhibitory", "MSN_1", "MSN_2", "MSN_3", "WM"
)

safe_colorblind_palette <- c(
  "#E7298A", "#A6761D", "#D95F02", "#E6AB02",
  "#66A61E", "#1B9E77", "#7570B3", "#666666"
)

names(safe_colorblind_palette) <- domain_levels

# extract metadata
df <- as.data.frame(colData(spe)) |>
  dplyr::select(sample_id, spatial_domains) |>
  mutate(
    sample_id = factor(sample_id, levels = plot_sample_order),
    spatial_domains = factor(spatial_domains, levels = domain_levels)
  )

# compute proportions
prop_df <- df |>
  count(sample_id, spatial_domains, name = "n_spots") |>
  group_by(sample_id) |>
  mutate(prop_spots = n_spots / sum(n_spots)) |>
  ungroup() |>
  complete(
    sample_id,
    spatial_domains,
    fill = list(n_spots = 0, prop_spots = 0)
  )

# faceted plot
p <- ggplot(prop_df, aes(x = sample_id, y = prop_spots, fill = spatial_domains)) +
  geom_col(width = 0.85) +
  facet_wrap(~ spatial_domains, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = safe_colorblind_palette, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Samples ordered by A-P position",
    y = "Proportion of spots",
    fill = "Spatial domain"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

plot_dir <- here("plots", "25_revisions", "07-prop_SpD")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(plot_dir, "proportion_SpD_AP_axis.pdf"), width = 6, height = 6)
print(p)
dev.off()