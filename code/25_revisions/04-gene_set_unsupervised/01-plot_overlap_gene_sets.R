####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
library(sessioninfo)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(spatialNAcUtils)
library(jaffelab)
library(projectR)
library(scater)
library(scran)
library(dittoSeq)
library(escheR)
library(getopt)
library(Seurat)

# Read human genes
x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', "human_NAc", paste0('nmf_results.rds')))
human_genes <- rownames(x@w)
rm(x)

# Read case control cocaine genes
x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', "rat_case_control_acute", paste0('nmf_results.rds')))
cocaine_acute_genes <- rownames(x@w)
rm(x)

x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', "rat_case_control_repeated", paste0('nmf_results.rds')))
cocaine_repeated_genes <- rownames(x@w)
rm(x)

# Read case control cocaine genes
x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', "rat_case_control_morphine_acute", paste0('nmf_results.rds')))
morphine_acute_genes <- rownames(x@w)
rm(x)

x <- readRDS(file=here::here('processed-data','16_transfer_learning','01_process_reference', 'RCppML', "rat_case_control_morphine_repeated", paste0('nmf_results.rds')))
morphine_repeated_genes <- rownames(x@w)
rm(x)

# Plot directory
plotDir <- here::here("plots", "25_revisions", "04-gene_set_unsupervised")
## ----------------------------
## 1) Bar plot of total genes
## ----------------------------
gene_counts <- tibble(
  dataset = c(
    "Human NAc",
    "Cocaine acute",
    "Cocaine repeated",
    "Morphine acute",
    "Morphine repeated"
  ),
  n_genes = c(
    length(human_genes),
    length(cocaine_acute_genes),
    length(cocaine_repeated_genes),
    length(morphine_acute_genes),
    length(morphine_repeated_genes)
  )
) %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c(
        "Human NAc",
        "Cocaine acute",
        "Cocaine repeated",
        "Morphine acute",
        "Morphine repeated"
      )
    )
  )

p_bar <- ggplot(gene_counts, aes(x = dataset, y = n_genes)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(n_genes, big.mark = ",")), vjust = -0.4, size = 4) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Number of genes",
    title = "Input genes used for factorization"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  expand_limits(y = max(gene_counts$n_genes) * 1.08)

pdf(file.path(plotDir, "gene_counts_barplot.pdf"), width = 5, height = 3.5)
print(p_bar + theme(legend.position = "none"))
dev.off()

library(ComplexHeatmap)

nonhuman_gene_sets <- list(
  cocaine_acute = unique(cocaine_acute_genes),
  cocaine_repeated = unique(cocaine_repeated_genes),
  morphine_acute = unique(morphine_acute_genes),
  morphine_repeated = unique(morphine_repeated_genes)
)

upset_input <- make_comb_mat(nonhuman_gene_sets)
upset_input <- upset_input[comb_size(upset_input) >= 30]



pdf(file.path(plotDir, "nonhuman_gene_sets_upset.pdf"), width = 6, height = 3)
draw(
  UpSet(
    upset_input,
    set_order = c(
      "cocaine_acute",
      "cocaine_repeated",
      "morphine_acute",
      "morphine_repeated"
    ),
    comb_order = order(comb_size(upset_input), decreasing = TRUE),
    top_annotation = upset_top_annotation(
      upset_input,
      add_numbers = TRUE, height = unit(3.5, "cm"), numbers_rot = 0, width = unit(5, "cm")
    ),
    right_annotation = upset_right_annotation(upset_input)
  )
)
dev.off()