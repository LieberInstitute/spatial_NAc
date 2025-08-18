library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(scater)
library(reshape2)
library(here)

set.seed(123)

opt <- list()
opt$data <- "human_NAc"

# Load the SPE data and add cluster information
nmf_dir <- here::here('processed-data', '16_transfer_learning')
spe <- readRDS(file.path(nmf_dir, "02_target_projections", opt$data, paste0("spe_NMF.rds")))
clusters_resFile <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["domain"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_resFile), by = 'key') |>
    pull(cluster) |>
    as.factor()

# Specify the sample order from anterior to posterior 
sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522", "Br8492", "Br8325", "Br8667", "Br3942")

plot_df <- colData(spe)
select_NMFs <- c("nmf38", "nmf10", "nmf3", "nmf7", "nmf39", "nmf4", "nmf25")
plot_df <- plot_df[ ,c(c("sample_id", "domain"), select_NMFs)]
plot_df$domain <- as.character(plot_df$domain)
plot_df <- plot_df[plot_df$domain %in% c("MSN 1", "MSN 2", "MSN 3"), ]
plot_df$domain <- factor(plot_df$domain, levels = c("MSN 1", "MSN 2", "MSN 3"))
plot_df$sample_id <- as.character(plot_df$sample_id)
plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)

library(tidyverse)
plot_df_long <- data.frame(plot_df) %>%
  pivot_longer(
    cols = starts_with("nmf"), 
    names_to = "nmf_factor", 
    values_to = "value"
  )
plotDir <- here::here('plots', '16_transfer_learning', '02_target_projections', opt$data)

pdf(file.path(plotDir, "MSN_domain_factors_along_AP_axis.pdf"), width = 10, height = 3)
# Loop over NMF factors
for (factor in unique(plot_df_long$nmf_factor)) {
  p <- ggplot(filter(plot_df_long, nmf_factor == factor), 
              aes(x = sample_id, y = value, fill = domain)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(
      values = c(
        "MSN 1" = "#66A61E", 
        "MSN 2" = "#1B9E77", 
        "MSN 3" = "#7570B3"
      )
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      strip.text = element_text(size = 10)
    ) +
    labs(
      title = paste("Grouped Boxplot for", factor),
      fill = "Domain"
    )
  
  print(p)  # Draw the plot on the PDF page
}
dev.off()