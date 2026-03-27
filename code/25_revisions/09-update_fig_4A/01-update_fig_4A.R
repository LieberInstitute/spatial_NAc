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
library("ggpubr")
library("rstatix")

# Specify the output directory
plotDir <- here("plots", "25_revisions", "09-update_fig_4A")

# Read in SPE object
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")
spe <- loadHDF5SummarizedExperiment(spe_dir_out)
spe <- spe[, !is.na(spe[["spatial_domains"]])]

safe_colorblind_palette <- c(
  "#66A61E", "#1B9E77", "#7570B3", "#E7298A",
  "#D95F02", "#E6AB02", "#666666", "#A6761D"
)

spe$spatial_domains <- factor(
  spe$spatial_domains,
  levels = c(
    "MSN_1", "MSN_2", "MSN_3", "D1 islands",
    "Excitatory", "Inhibitory", "WM", "Endothelial/Ependymal"
  )
)

# -----------------------------
# Build long data frame
# -----------------------------
plot_df <- data.frame(colData(spe))
plot_df <- plot_df[, c(
  "sample_id", "slide_num", "array_num", "spatial_domains", "exclude_overlapping",
  "DRD1_MSN_A", "DRD2_MSN_A", "DRD2_MSN_B", "DRD1_MSN_C", "DRD1_MSN_B", "DRD1_MSN_D",
  "Excitatory", "Inh_A", "Inh_B", "Inh_C", "Inh_D", "Inh_E", "Inh_F",
  "Astrocyte_A", "Astrocyte_B", "Oligo", "OPC", "Microglia", "Endothelial", "Ependymal"
)]

plot_df$capture_area <- paste(plot_df$slide_num, plot_df$array_num, sep = "_")

# optional
plot_df <- plot_df %>%
  filter(!is.na(exclude_overlapping)) %>%
  filter(!exclude_overlapping)

plot_df <- reshape2::melt(
  plot_df,
  id.vars = c("sample_id", "slide_num", "array_num", "capture_area", "spatial_domains"),
  variable.name = "cell_type",
  value.name = "RCTD_weights"
)

colnames(plot_df)[colnames(plot_df) == "spatial_domains"] <- "spatial_domain"

plot_df <- plot_df[plot_df$spatial_domain %in% c("MSN_1", "MSN_2", "MSN_3"), ]
plot_df$cell_type <- gsub("_", " ", plot_df$cell_type)

select_cell_types <- c(
  "DRD1 MSN A", "DRD1 MSN B", "DRD1 MSN C", "DRD1 MSN D",
  "DRD2 MSN A", "DRD2 MSN B", "Oligo", "Astrocyte A"
)

plot_df_subset <- plot_df[plot_df$cell_type %in% select_cell_types, ]

plot_df_subset$cell_type <- factor(plot_df_subset$cell_type, levels = select_cell_types)
plot_df_subset$spatial_domain <- factor(
  plot_df_subset$spatial_domain,
  levels = c("MSN_1", "MSN_2", "MSN_3")
)

# -----------------------------
# Aggregate by capture area for testing only
# -----------------------------
test_df <- plot_df_subset %>%
  group_by(sample_id, capture_area, cell_type, spatial_domain) %>%
  summarise(
    RCTD_weights = mean(RCTD_weights, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Pairwise Wilcoxon tests within each cell type
# -----------------------------
stat_test <- test_df %>%
  group_by(cell_type) %>%
  wilcox_test(RCTD_weights ~ spatial_domain,
              comparisons = list(
                c("MSN_1", "MSN_2"),
                c("MSN_1", "MSN_3"),
                c("MSN_2", "MSN_3")
              )) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# keep only significant comparisons
stat_test <- stat_test %>%
  filter(p.adj < 0.05)

# -----------------------------
# Add y positions for the plotted data
# -----------------------------
# This uses the spot-level data to place brackets above the boxplots
stat_test <- stat_test %>%
  add_xy_position(
    x = "cell_type",
    group = "spatial_domain",
    dodge = 0.8,
    data = plot_df_subset
  )

# -----------------------------
# Plot
# -----------------------------
pdf(
  file.path(plotDir, "RCTD_weights_boxplot_marker_genes_select_cell_types_sigbars.pdf"),
  width = 5.2, height = 7.2
)

p <- ggplot(
  plot_df_subset,
  aes(x = cell_type, y = RCTD_weights)
) +
  geom_boxplot(
    aes(fill = spatial_domain),
    outlier.shape = NA,
    width = 0.55,
    position = position_dodge(width = 0.8)
  ) +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.01,
    coord.flip = TRUE,
    hide.ns = TRUE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("#66A61E", "#1B9E77", "#7570B3")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  guides(fill = guide_legend(title = "Spatial Domain")) +
  theme_pubr() +
  coord_flip(ylim = c(0, 1.42), clip = "off") +
  xlab("") +
  ylab("RCTD weights") +
  labs(caption = "* adj. P < 0.05; ** adj. P < 0.01; *** adj. P < 0.001")
print(p)
dev.off()