rm(list = ls())

library(liana)
library(tidyverse)
library(here)

# -----------------------------
# Input / output
# -----------------------------
res_dir <- here("processed-data", "22_gene_risk_LR_analysis", "01_liana")
plot_dir <- here("plots", "25_revisions", "12-non_neuronal_MSN", "LR_glia_to_MSN")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

liana_res <- readRDS(file.path(res_dir, "liana_consensus.rds"))

# keep only consensus interactions
liana_trunc <- liana_res %>%
  filter(aggregate_rank <= 0.01)

# receiver populations of interest
receiver_levels <- c("DRD2_MSN_A", "DRD1_MSN_A", "DRD1_MSN_C")

# -----------------------------
# Helper function
# -----------------------------
plot_sender_to_receivers <- function(
    liana_df,
    sender_cell,
    receiver_cells = receiver_levels,
    top_n = 8,
    outdir = plot_dir
) {
  
  df_sub <- liana_df %>%
    filter(
      source == sender_cell,
      target %in% receiver_cells
    ) %>%
    mutate(
      LR_pair = paste0(ligand.complex, " -> ", receptor.complex)
    )
  
  if (nrow(df_sub) == 0) {
    stop("No interactions found for sender: ", sender_cell)
  }
  
  # pick top interactions separately for each receiver
  top_hits <- df_sub %>%
    group_by(target) %>%
    arrange(aggregate_rank, desc(sca.LRscore), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # keep the union of selected LR pairs
  selected_pairs <- unique(top_hits$LR_pair)
  
  df_plot <- df_sub %>%
    filter(LR_pair %in% selected_pairs)
  
  # order rows by best aggregate rank across displayed interactions
  row_order <- df_plot %>%
    group_by(LR_pair) %>%
    summarise(
      best_rank = min(aggregate_rank, na.rm = TRUE),
      best_score = max(sca.LRscore, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_rank, desc(best_score)) %>%
    pull(LR_pair)
  
  df_plot <- df_plot %>%
    mutate(
      target = factor(target, levels = receiver_cells),
      LR_pair = factor(LR_pair, levels = rev(row_order))
    )
  
  # save selected top hits
  write.csv(
    top_hits %>%
      select(
        source, target, ligand.complex, receptor.complex, LR_pair,
        aggregate_rank, natmi.edge_specificity, sca.LRscore,
        connectome.weight_sc, logfc.logfc_comb, cellphonedb.pvalue
      ) %>%
      arrange(target, aggregate_rank),
    file.path(outdir, paste0("top", top_n, "_", sender_cell, "_to_receivers.csv")),
    row.names = FALSE
  )
  
  # save displayed plot table
  write.csv(
    df_plot %>%
      select(
        source, target, ligand.complex, receptor.complex, LR_pair,
        aggregate_rank, natmi.edge_specificity, sca.LRscore,
        connectome.weight_sc, logfc.logfc_comb, cellphonedb.pvalue
      ) %>%
      arrange(LR_pair, target),
    file.path(outdir, paste0("displayed_", sender_cell, "_to_receivers.csv")),
    row.names = FALSE
  )
  
  p <- ggplot(
    df_plot,
    aes(
      x = target,
      y = LR_pair,
      size = natmi.edge_specificity,
      color = sca.LRscore
    )
  ) +
    geom_point() +
    scale_color_viridis_c(
      option = "viridis",
      name = "Expression\nmagnitude"
    ) +
    scale_size(
      range = c(2, 6),
      name = "Interaction\nspecificity"
    ) +
    labs(
      x = "Receiver",
      y = NULL,
      title = paste0(sender_cell, " -> MSN interactions")
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  plot_height <- max(5, 0.28 * length(unique(df_plot$LR_pair)) + 1.5)
  
  ggsave(
    filename = file.path(outdir, paste0("dotplot_", sender_cell, "_to_receivers.pdf")),
    plot = p,
    width = 5,
    height = plot_height
  )
  
  return(list(
    plot = p,
    top_hits = top_hits,
    displayed = df_plot
  ))
}

# -----------------------------
# Run plots
# -----------------------------
astro_res <- plot_sender_to_receivers(
  liana_df = liana_trunc,
  sender_cell = "Astrocyte_A",
  receiver_cells = receiver_levels,
  top_n = 15,
  outdir = plot_dir
)

oligo_res <- plot_sender_to_receivers(
  liana_df = liana_trunc,
  sender_cell = "Oligo",
  receiver_cells = receiver_levels,
  top_n = 15,
  outdir = plot_dir
)