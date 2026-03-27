#### code for plotting Reactome GSEA results for selected MCP factors

library(dplyr)
library(ggplot2)
library(reactome.db)
library(org.Hs.eg.db)
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(gridExtra)

set.seed(123)

# ----------------------------
# I/O
# ----------------------------
dat_dir  <- here("code", "06_deploy_app", "spe_shiny")
plot_dir <- here::here("plots", "25_revisions", "01_MCP_GSEA")
save_dir <- here::here("processed-data", "25_revisions", "01_MCP_GSEA")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load data
# ----------------------------
spe <- loadHDF5SummarizedExperiment(dat_dir)

# ----------------------------
# Prepare expression matrix
# ----------------------------
expr_mat <- logcounts(spe)
geneData <- rowData(spe)

rownames(expr_mat) <- geneData$gene_id

# retain unique protein-coding non-mito/non-ribosomal genes
gene_keep <- !duplicated(geneData$gene_name) &
  geneData$gene_type == "protein_coding" &
  !grepl("^MT-", geneData$gene_name) &
  !grepl("^RP[SL]", geneData$gene_name)

geneData <- geneData[gene_keep, ]
expr_mat <- expr_mat[rownames(expr_mat) %in% geneData$gene_id, , drop = FALSE]

# ----------------------------
# Compute gene x MCP correlations
# ----------------------------
pattern_cols <- grep("^meringue_cluster_", colnames(colData(spe)), value = TRUE)

compute_cor_matrix <- function(expr, patterns, spe_obj) {
  cor_mat <- matrix(NA_real_, nrow = nrow(expr), ncol = length(patterns))
  rownames(cor_mat) <- rownames(expr)
  colnames(cor_mat) <- patterns

  for (pat in patterns) {
    message("Computing correlations for ", pat)
    pat_values <- colData(spe_obj)[[pat]]
    non_na_idx <- which(!is.na(pat_values))

    expr_sub <- expr[, non_na_idx, drop = FALSE]
    pat_sub  <- pat_values[non_na_idx]

    expr_sub <- as.matrix(expr_sub)

    cors <- apply(expr_sub, 1, function(x) cor(x, pat_sub, method = "pearson"))
    cor_mat[, pat] <- cors
  }

  cor_mat
}

cor_mat <- compute_cor_matrix(expr_mat, pattern_cols, spe)
cor_mat <- cor_mat[rowSums(is.na(cor_mat)) == 0, , drop = FALSE]

# ----------------------------
# Prep Reactome
# ----------------------------
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo", xx)]

reactome_ids <- as.list(reactomePATHID2EXTID)
reactome.h <- reactome.h[intersect(names(reactome.h), names(reactome_ids))]
x.h <- reactome_ids[names(reactome.h)]

reactome.h <- gsub("Homo sapiens: ", "", reactome.h)
names(x.h) <- reactome.h

# ----------------------------
# Helper to read one MCP result file
# ----------------------------
read_mcp_results <- function(mcp_short, save_dir) {
  f <- file.path(save_dir, paste0(mcp_short, "_reactome_results.csv"))
  if (!file.exists(f)) {
    stop("Could not find results file: ", f)
  }
  read.csv(f, row.names = 1, stringsAsFactors = FALSE)
}

# ----------------------------
# Custom plotting function
# ----------------------------
gseaPlotMCP <- function(mcp_name, mcp_short, mcp_terms, cor_mat, x.h, save_dir) {
  results <- read_mcp_results(mcp_short, save_dir)

  tmp <- results[results$pathway %in% mcp_terms, , drop = FALSE]
  if (nrow(tmp) == 0) {
    stop("None of the requested pathways were found for ", mcp_short)
  }

  # preserve requested order
  tmp$pathway <- factor(tmp$pathway, levels = mcp_terms)
  tmp <- tmp[order(tmp$pathway), , drop = FALSE]
  tmp$pathway <- as.character(tmp$pathway)

  # positively correlated genes only, matching enrichment script
  non0.mcp <- rownames(cor_mat)[cor_mat[, mcp_name] > 0]

  non0.id <- mapIds(
    org.Hs.eg.db,
    keys = non0.mcp,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  names(non0.id) <- non0.mcp
  non0.id <- non0.id[!is.na(non0.id)]

  mcp.stats <- cor_mat[names(non0.id), mcp_name]
  names(mcp.stats) <- non0.id
  mcp.stats <- sort(mcp.stats, decreasing = TRUE)

  plot_df <- do.call(rbind, lapply(seq_len(nrow(tmp)), function(i) {
    term1 <- tmp$pathway[i]
    data.frame(
      id = names(mcp.stats),
      weight = as.numeric(mcp.stats),
      rank = seq_along(mcp.stats),
      term = term1,
      present = names(mcp.stats) %in% x.h[[term1]],
      y = i - 1,
      yend = i,
      stringsAsFactors = FALSE
    )
  }))

  y_labels <- sapply(seq_len(nrow(tmp)), function(i) {
    paste0(
      "NES= ", round(tmp$NES[i], 2),
      "\n(", tmp$size[i], "/", length(x.h[[tmp$pathway[i]]]), ")",
      "\npadj= ", format(tmp$padj[i], scientific = TRUE, digits = 2)
    )
  })

  ggplot(dplyr::filter(plot_df, present), aes(x = rank, xend = rank, y = y, yend = yend, color = term)) +
    geom_segment(linewidth = 0.6) +
    scale_x_continuous(
      limits = c(0, max(plot_df$rank)),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(0.5, nrow(tmp) - 0.5, by = 1),
      labels = y_labels,
      expand = c(0, 0)
    ) +
    theme_minimal() +
    ggtitle(mcp_short) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
}

# ----------------------------
# Final pathway sets
# ----------------------------

mcp1.terms <- c(
  "Transmission across Chemical Synapses",
  "Signaling by GPCR",
  "DARPP-32 events",
  "Opioid Signalling",
  "GABA synthesis, release, reuptake and degradation",
  "Post NMDA receptor activation events",
  "EPH-Ephrin signaling"
)

mcp3.terms <- c(
  "Transmission across Chemical Synapses",
  "Signaling by GPCR",
  "Post NMDA receptor activation events",
  "Glutamate binding, activation of AMPA receptors and synaptic plasticity",
  "Opioid Signalling",
  "L1CAM interactions",
  "Neurotransmitter release cycle"
)

mcp4.terms <- c(
  "Signaling by GPCR",
  "G alpha (i) signalling events",
  "PLC beta mediated events",
  "DAG and IP3 signaling",
  "Ca-dependent events",
  "Transmission across Chemical Synapses",
  "Opioid Signalling"
)

# ----------------------------
# Make plots
# ----------------------------
p1 <- gseaPlotMCP(
  mcp_name = "meringue_cluster_1",
  mcp_short = "mcp1",
  mcp_terms = mcp1.terms,
  cor_mat = cor_mat,
  x.h = x.h,
  save_dir = save_dir
)

p3 <- gseaPlotMCP(
  mcp_name = "meringue_cluster_3",
  mcp_short = "mcp3",
  mcp_terms = mcp3.terms,
  cor_mat = cor_mat,
  x.h = x.h,
  save_dir = save_dir
)

p4 <- gseaPlotMCP(
  mcp_name = "meringue_cluster_4",
  mcp_short = "mcp4",
  mcp_terms = mcp4.terms,
  cor_mat = cor_mat,
  x.h = x.h,
  save_dir = save_dir
)

# ----------------------------
# Save combined figure
# ----------------------------
ggsave(
  filename = file.path(plot_dir, "MCP_custom_gsea_plots_selected.pdf"),
  plot = gridExtra::grid.arrange(p1, p3, p4, ncol = 1),
  bg = "white",
  height = 18,
  width = 12,
  units = "in"
)

# ----------------------------
# Save individual panels too
# ----------------------------
ggsave(
  filename = file.path(plot_dir, "MCP1_custom_gsea_plot.pdf"),
  plot = p1,
  bg = "white",
  height = 6,
  width = 12,
  units = "in"
)

ggsave(
  filename = file.path(plot_dir, "MCP3_custom_gsea_plot.pdf"),
  plot = p3,
  bg = "white",
  height = 6,
  width = 12,
  units = "in"
)

ggsave(
  filename = file.path(plot_dir, "MCP4_custom_gsea_plot.pdf"),
  plot = p4,
  bg = "white",
  height = 6,
  width = 12,
  units = "in"
)