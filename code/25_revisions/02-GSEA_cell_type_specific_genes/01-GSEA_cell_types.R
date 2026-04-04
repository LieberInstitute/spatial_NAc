#### Reactome pathway enrichment for upregulated cell type DEGs using gprofiler2
#### Save minimal output table

library(dplyr)
library(gprofiler2)
library(here)

set.seed(123)

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
dat_dir <- here("processed-data", "10_post_clustering_analysis", "02_spatial_registration_sn")
sn_registration <- readRDS(file.path(dat_dir, "sn_cellType_registration.rds"))
enrichment <- sn_registration$enrichment

# ------------------------------------------------------------------
# Define cell types from enrichment table
# ------------------------------------------------------------------
tstat_cols <- grep("^t_stat_", colnames(enrichment), value = TRUE)
cell_types <- sub("^t_stat_", "", tstat_cols)

# ------------------------------------------------------------------
# Define background = all detected genes in the study
# ------------------------------------------------------------------
background_genes <- enrichment$gene
background_genes <- background_genes[!is.na(background_genes) & background_genes != ""]
background_genes <- unique(background_genes)

# ------------------------------------------------------------------
# Function to run Reactome enrichment for one cell type
# Query genes: FDR < 0.05 and logFC > 0.5
# Multiple testing correction: BH/FDR via g:Profiler
# ------------------------------------------------------------------
run_gprofiler_reactome <- function(enrichment_df,
                                   cell_type,
                                   fdr_cutoff = 0.05,
                                   logfc_cutoff = 0.5) {
  
  fdr_col <- paste0("fdr_", cell_type)
  logfc_col <- paste0("logFC_", cell_type)
  
  query_df <- enrichment_df %>%
    dplyr::select(gene, !!fdr_col, !!logfc_col) %>%
    dplyr::rename(
      fdr = !!fdr_col,
      logFC = !!logfc_col
    ) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::filter(!is.na(fdr), !is.na(logFC)) %>%
    dplyr::filter(fdr < fdr_cutoff, logFC > logfc_cutoff)
  
  query_genes <- unique(query_df$gene)
  
  if (length(query_genes) < 2) {
    message("Skipping ", cell_type, ": fewer than 2 genes pass thresholds")
    return(NULL)
  }
  
  gp_res <- gost(
    query = query_genes,
    organism = "hsapiens",
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = FALSE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = FALSE,
    user_threshold = 0.05,
    correction_method = "fdr",
    custom_bg = background_genes,
    sources = c("REAC")
  )
  
  if (is.null(gp_res) || is.null(gp_res$result) || nrow(gp_res$result) == 0) {
    message("No Reactome enrichment for ", cell_type)
    return(NULL)
  }
  
  out <- gp_res$result %>%
    dplyr::transmute(
      cell_type = cell_type,
      pathway = term_name,
      p_adj_BH = p_value,
      intersection_size = intersection_size,
      term_size = term_size
    ) %>%
    dplyr::filter(p_adj_BH < 0.05) %>%
    dplyr::arrange(p_adj_BH)
  
  return(out)
}

# ------------------------------------------------------------------
# Run across all cell types
# ------------------------------------------------------------------
reactome_list <- lapply(cell_types, function(ct) {
  message("Running Reactome enrichment for ", ct)
  run_gprofiler_reactome(enrichment, ct)
})

# Remove NULL results
reactome_list <- reactome_list[!vapply(reactome_list, is.null, logical(1))]

# Bind results
if (length(reactome_list) > 0) {
  reactome_df <- bind_rows(reactome_list) %>%
    dplyr::arrange(cell_type, p_adj_BH)
} else {
  reactome_df <- data.frame(
    cell_type = character(),
    pathway = character(),
    p_adj_BH = numeric(),
    intersection_size = integer(),
    term_size = integer()
  )
}

# ------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------
out_file <- file.path(dat_dir, "gprofiler_REACTOME_celltypes_BH_minimal.csv")
write.csv(reactome_df, out_file, row.names = FALSE)

# ------------------------------------------------------------------
# Optional preview
# ------------------------------------------------------------------
print(head(reactome_df))
cat("Saved results to:", out_file, "\n")