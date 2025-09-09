library("SpatialExperiment")
library("HDF5Array")
library("lobstr")
library("here")
library("sessioninfo")
library("Matrix")
library("dplyr")

spe_dir_in <- here(
    "processed-data", "05_harmony_BayesSpace", "04-preprocess_and_harmony", "spe_harmony.rds"
)
spe_dir_out <- here("code", "06_deploy_app", "spe_shiny")

spe <- readRDS(spe_dir_in)
spe

#class: SpatialExperiment 
#dim: 22375 176013 
#metadata(0):
#assays(3): counts logcounts binomial_deviance_residuals
#rownames(22375): AL627309.1 AL627309.5 ... AC007325.4 AC007325.2
#rowData names(15): source type ... per.block binomial_deviance
#colnames(176013): AAACAACGAATAGTTC-1_V11D01-384_A1
#  AAACAAGTATCTCCCA-1_V11D01-384_A1 ... TTGTTTGTATTACACG-1_V13M13-362_C1
#  TTGTTTGTGTAAATTC-1_V13M13-362_C1
#colData names(93): sample_id in_tissue ... row col
#reducedDimNames(22): 10x_pca 10x_tsne ... TSNE_perplexity80.HARMONY
#  UMAP.HARMONY
#mainExpName: NULL
#altExpNames(0):
#spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#imgData names(4): sample_id image_id data scaleFactor

#   All images should be low-res (if others existed, we'd drop them)
stopifnot(all(imgData(spe)$image_id == "lowres"))

#   Only keep logcounts
assays(spe) <- list(logcounts = assays(spe)$logcounts)

#############################Add new metadata######################################
remove_cols <- grepl("SNN_k10", colnames(colData(spe)))
colData(spe) <- colData(spe)[ ,!remove_cols]

spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)
sample_ids <- levels(spe$sample_id)

opt <- list()
opt$marker_genes <- TRUE
RCTD_list <- lapply(sample_ids, function(isample){
    RCTD_dir <- here::here("processed-data", "08_spot_deconvo", "01_RCTD", isample)
    if(opt$marker_genes){
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD_markers.rds"))
    }else{
        myRCTD <- readRDS(file.path(RCTD_dir, "results_RCTD.rds"))
    }
    myRCTD
})

weights <- lapply(RCTD_list, function(myRCTD){
    my_results = myRCTD@results
    my_weights = lapply(my_results, function(x) x$all_weights)
    my_weights_df <- data.frame(do.call(rbind, my_weights))
    my_weights_df
})
weights <- data.frame(do.call(rbind, weights))

coords <- lapply(RCTD_list, function(myRCTD){
    my_coords = myRCTD@spatialRNA@coords
    my_coords
})
coords <- data.frame(do.call(rbind, coords))

spe <- spe[ ,colnames(spe) %in% rownames(coords)]
spe <- spe[ ,match(rownames(coords), colnames(spe))]
colData(spe) <- cbind(colData(spe), weights)

# Add the unmerged clusters
for(K in c(3:15)){
    unmerged_res_path <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast",
     "random_start_3", paste0("PRECAST_k",K ,".csv"))
    precast_results <- unmerged_res_path |>
        read.csv() |>
        as_tibble() |>
        select(c(key, cluster))
    colnames(precast_results) <- c("key", paste0("PRECAST_K",K ,"_clusters"))
    temp <- colnames(spe)
    colData(spe) <- colData(spe) |>
        as_tibble() |>
        left_join(precast_results, by = "key") |>
        DataFrame()
    colnames(spe) <- temp
}

# Add the final spatial domains
clusters_file <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters", "precast_clusters.csv")
spe[["spatial_domains"]] = colData(spe) |>
    as_tibble() |>
    left_join(read.csv(clusters_file), by = 'key') |>
    pull(cluster) |>
    as.factor()

spe <- spe[ ,!is.na(spe[["spatial_domains"]])]
# Add MERINGUE patters
meringue_res_path <- here("processed-data", "14_MSN_factorization", "01_meringue", "meringue_consensus_patterns.rds")
meringue_df <- readRDS(meringue_res_path)
key_spe  <- rownames(colData(spe))
key_mer  <- rownames(meringue_df)
idx <- match(key_spe, key_mer)
n_tot <- length(idx)
n_hit <- sum(!is.na(idx))
message(sprintf("Matched %d / %d spots (%.1f%%); %d unmatched get NA.",
                n_hit, n_tot, 100*n_hit/n_tot, n_tot - n_hit))
aligned <- meringue_df[idx, , drop = FALSE] # adds NA rows wherever idx is NA
rownames(aligned) <- rownames(colData(spe)) # keep spe order
colData(spe) <- cbind(colData(spe), aligned)

# Add human NMF factors
nmf_dir <- here::here('processed-data', '16_transfer_learning')
spe_human_NMF <- readRDS(file.path(nmf_dir, "02_target_projections/human_NAc/spe_NMF.rds"))
nmf_df <- colData(spe_human_NMF)
nmf_df <- nmf_df[ ,grepl("nmf", colnames(nmf_df))]
colnames(nmf_df) <- paste0("human_snRNA_", colnames(nmf_df))
key_nmf <- rownames(nmf_df)
idx <- match(key_spe, key_nmf)
aligned_nmf <- nmf_df[idx, , drop = FALSE] # adds NA rows wherever idx is NA
rownames(aligned_nmf) <- rownames(colData(spe)) # keep spe order
colData(spe) <- cbind(colData(spe), aligned_nmf)
rm(spe_human_NMF)
# Add rat case control volitional morphine analysis
spe_morphine_repeated_NMF <- readRDS(file.path(nmf_dir, "04_cellType_NMF/rat_case_control_morphine_repeated/spe_NMF_30.rds"))
morphine_repeated_df <- colData(spe_morphine_repeated_NMF)
rm(spe_morphine_repeated_NMF)
spe_morphine_acute_NMF <- readRDS(file.path(nmf_dir, "04_cellType_NMF/rat_case_control_morphine_acute/spe_NMF_30.rds"))
morphine_acute_df <- colData(spe_morphine_acute_NMF)
rm(spe_morphine_acute_NMF)
spe_cocaine_acute_NMF <- readRDS(file.path(nmf_dir, "04_cellType_NMF/rat_case_control_cocaine_acute/spe_NMF_30.rds"))
cocaine_acute_df <- colData(spe_cocaine_acute_NMF)
rm(spe_cocaine_acute_NMF)

morphine_repeated_df <- morphine_repeated_df[ ,grepl("nmf", colnames(morphine_repeated_df))]
morphine_acute_df <- morphine_acute_df[ ,grepl("nmf", colnames(morphine_acute_df))]
cocaine_acute_df <- cocaine_acute_df[ ,grepl("nmf", colnames(cocaine_acute_df))]

colnames(morphine_repeated_df) <- paste0("morphine_volitional_", colnames(morphine_repeated_df))
colnames(morphine_acute_df) <- paste0("morphine_acute_", colnames(morphine_acute_df))
colnames(cocaine_acute_df) <- paste0("cocaine_acute_", colnames(cocaine_acute_df))

drug_nmf_df <- cbind(morphine_repeated_df, morphine_acute_df, cocaine_acute_df)
key_drug <- rownames(drug_nmf_df)
idx <- match(key_spe, key_drug)
aligned_drug <- drug_nmf_df[idx, , drop = FALSE] # adds NA rows wherever idx is NA
rownames(aligned_drug) <- rownames(colData(spe)) # keep spe order
colData(spe) <- cbind(colData(spe), aligned_drug)


###################################################################################
# Make some plots to check this object
plotDir <- here("plots", "06_deploy_app")

library(spatialLIBD)
plot_sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522", "Br8492", "Br8325", "Br8667", "Br3942")

spe[["spatial_domains"]] <- as.character(spe[["spatial_domains"]])
spe[["spatial_domains"]][spe[["spatial_domains"]] == "MSN 1"] <- "MSN_1"
spe[["spatial_domains"]][spe[["spatial_domains"]] == "MSN 2"] <- "MSN_2"
spe[["spatial_domains"]][spe[["spatial_domains"]] == "MSN 3"] <- "MSN_3"

levels(spe[["spatial_domains"]]) <- c("MSN_1", "MSN_2", "MSN_3", "D1 islands", "Inhibitory", "Excitatory", "Endothelial/Ependymal", "WM")
vis_grid_clus(spe, clustervar = "spatial_domains", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "spatial_domains.pdf"), sort_clust = FALSE)

vis_grid_gene(spe, geneid = "meringue_cluster_1", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "meringue_cluster_1.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_2", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "meringue_cluster_2.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_3", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "meringue_cluster_3.pdf"))

vis_grid_gene(spe, geneid = "meringue_cluster_4", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "meringue_cluster_4.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf3", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf3.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf4", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf4.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf7", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf7.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf10", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf10.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf39", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf39.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf34", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf34.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf35", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf35.pdf"))

vis_grid_gene(spe, geneid = "human_snRNA_nmf44", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "human_snRNA_nmf44.pdf"))

vis_grid_gene(spe, geneid = "morphine_volitional_nmf12", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "morphine_volitional_nmf12.pdf"))

vis_grid_gene(spe, geneid = "morphine_volitional_nmf3", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "morphine_volitional_nmf3.pdf"))

vis_grid_gene(spe, geneid = "DRD1_MSN_B", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "DRD1_MSN_B.pdf"))

vis_grid_gene(spe, geneid = "DRD1_MSN_D", spatial = TRUE, sample_order = plot_sample_order, is_stitched = TRUE, 
pdf_file = file.path(plotDir, "DRD1_MSN_D.pdf"))


############################## Size considerations#################################
obj_size(spe)
# 978.00 MB
# (because logcounts are on disk, but Shiny appears to respect on-disk operations)

# (Comment out original gene-filtering steps, as the HDF5-backed SPE seems to
# behave well on Shiny, not loading all gene counts into memory)

# ## Check size in memory
# obj_size(spe)
# # 3.52 GB
# ## That's still too large for shinyapps

# ## Let's find out how many spots express each gene
# gene_expr <- rowSums(logcounts(spe) > 0)

# ## How many genes are expressed in 10 or more spots?
# table(gene_expr > 10)
# # FALSE  TRUE
# #  5561 25293

# ## 20 spots?
# table(gene_expr > 20)
# # FALSE  TRUE
# #  7053 23801

# ## 30 spots?
# table(gene_expr > 30)
# # FALSE  TRUE
# #  8018 22836

# ## 100 spots?
# table(gene_expr > 100)
# # FALSE  TRUE
# # 10794 20060

# ## 1000 spots?
# table(gene_expr > 1000)
# # FALSE  TRUE
# # 16018 14836

# ## 10000 spots?
# table(gene_expr > 10000)
# # FALSE  TRUE
# # 23085  7769

# obj_size(logcounts(spe))
# # 3.38 GB
# obj_size(logcounts(spe[gene_expr > 1000, ]))
# # 3.35 GB
# obj_size(logcounts(spe[gene_expr > 10000, ]))
# # 2.97 GB
# ## We have to go to genes expressed in at least 10,000 spots to drop under 3 GB
# ## for the logcounts

# ## Let's try first without doing so
# # spe <- spe[gene_expr > 10000, ]
# obj_size(spe)
# # 3.52 GB

# spe
# # class: SpatialExperiment
# # dim: 30854 135597
# # metadata(0):
# # assays(1): logcounts
# # rownames(30854): ENSG00000243485 ENSG00000238009 ... ENSG00000278817 ENSG00000277196
# # rowData names(7): source type ... gene_type gene_search
# # colnames(135597): AAACAACGAATAGTTC-1_V11D01-384_B1 AAACAAGTATCTCCCA-1_V11D01-384_B1 ...
# #   TTGTTTGTATTACACG-1_V13M06-378_D1 TTGTTTGTGTAAATTC-1_V13M06-378_D1
# # colData names(41): sample_id in_tissue ... scran_quick_cluster sizeFactor
# # reducedDimNames(3): 10x_pca 10x_tsne 10x_umap
# # mainExpName: NULL
# # altExpNames(0):
# # spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
# # imgData names(4): sample_id image_id data scaleFactor
saveHDF5SummarizedExperiment(spe, spe_dir_out, replace = TRUE)

## Reproducibility information
print("Reproducibility information:")
Sys.time() 
proc.time()

options(width = 120)
sessioninfo::session_info()

# [1] "Reproducibility information:"
# [1] "2025-09-02 15:07:16 EDT"
#      user    system   elapsed 
#  2155.004    62.427 12516.791 
# ─ Session info ──────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 (2023-10-31)
#  os       Rocky Linux 9.4 (Blue Onyx)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2025-09-02
#  pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc

# ─ Packages ──────────────────────────────────────────────────────────────────────
#  package                * version   date (UTC) lib source
#  abind                  * 1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
#  AnnotationDbi            1.64.1    2023-11-03 [1] Bioconductor
#  AnnotationHub            3.10.0    2023-10-24 [1] Bioconductor
#  attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.3.2)
#  beachmat                 2.18.0    2023-10-24 [1] Bioconductor
#  beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.3.2)
#  benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.3.2)
#  benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.3.2)
#  Biobase                * 2.62.0    2023-10-24 [1] Bioconductor
#  BiocFileCache            2.10.1    2023-10-26 [1] Bioconductor
#  BiocGenerics           * 0.48.1    2023-11-01 [1] Bioconductor
#  BiocIO                   1.12.0    2023-10-24 [1] Bioconductor
#  BiocManager              1.30.22   2023-08-08 [1] CRAN (R 4.3.2)
#  BiocNeighbors            1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
#  BiocParallel             1.36.0    2023-10-24 [1] Bioconductor
#  BiocSingular             1.18.0    2023-10-24 [1] Bioconductor
#  BiocVersion              3.18.1    2023-11-15 [1] Bioconductor
#  Biostrings               2.70.1    2023-10-25 [1] Bioconductor
#  bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
#  blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
#  bluster                  1.11.4    2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
#  bslib                    0.6.1     2023-11-28 [1] CRAN (R 4.3.2)
#  cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
#  cli                      3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
#  cluster                  2.1.6     2023-12-01 [1] CRAN (R 4.3.2)
#  codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
#  colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
#  config                   0.3.2     2023-08-30 [1] CRAN (R 4.3.2)
#  cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.3.2)
#  crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  curl                     5.2.0     2023-12-08 [1] CRAN (R 4.3.2)
#  data.table               1.14.10   2023-12-08 [1] CRAN (R 4.3.2)
#  DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
#  dbplyr                   2.4.0     2023-10-26 [1] CRAN (R 4.3.1)
#  DelayedArray           * 0.28.0    2023-10-24 [1] Bioconductor
#  DelayedMatrixStats       1.24.0    2023-10-24 [1] Bioconductor
#  digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
#  doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.3.2)
#  dotCall64                1.1-1     2023-11-28 [1] CRAN (R 4.3.2)
#  dplyr                  * 1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
#  dqrng                    0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
#  DT                       0.31      2023-12-09 [1] CRAN (R 4.3.2)
#  edgeR                    4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
#  ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
#  ExperimentHub            2.10.0    2023-10-24 [1] Bioconductor
#  fansi                    1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
#  farver                   2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
#  fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
#  fields                   15.2      2023-08-17 [1] CRAN (R 4.3.2)
#  filelock                 1.0.3     2023-12-11 [1] CRAN (R 4.3.2)
#  foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
#  generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb           * 1.38.1    2023-11-08 [1] Bioconductor
#  GenomeInfoDbData         1.2.11    2023-12-12 [1] Bioconductor
#  GenomicAlignments        1.38.0    2023-10-24 [1] Bioconductor
#  GenomicRanges          * 1.54.1    2023-10-29 [1] Bioconductor
#  ggbeeswarm               0.7.2     2023-04-29 [1] CRAN (R 4.3.2)
#  ggplot2                  3.5.1     2024-04-23 [1] CRAN (R 4.3.2)
#  ggrepel                  0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
#  glue                     1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
#  golem                    0.4.1     2023-06-05 [1] CRAN (R 4.3.2)
#  gridExtra                2.3       2017-09-09 [1] CRAN (R 4.3.2)
#  gtable                   0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
#  HDF5Array              * 1.30.0    2023-10-24 [1] Bioconductor
#  here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
#  htmltools                0.5.7     2023-11-03 [1] CRAN (R 4.3.2)
#  htmlwidgets              1.6.4     2023-12-06 [1] CRAN (R 4.3.2)
#  httpuv                   1.6.13    2023-12-06 [1] CRAN (R 4.3.2)
#  httr                     1.4.7     2023-08-15 [1] CRAN (R 4.3.1)
#  igraph                   2.0.3     2024-03-13 [1] CRAN (R 4.3.2)
#  interactiveDisplayBase   1.40.0    2023-10-24 [1] Bioconductor
#  IRanges                * 2.36.0    2023-10-24 [1] Bioconductor
#  irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
#  iterators                1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
#  jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.3.0)
#  jsonlite                 1.8.8     2023-12-04 [1] CRAN (R 4.3.2)
#  KEGGREST                 1.42.0    2023-10-24 [1] Bioconductor
#  labeling                 0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
#  later                    1.3.2     2023-12-06 [1] CRAN (R 4.3.2)
#  lattice                  0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
#  lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
#  lifecycle                1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
#  limma                    3.58.1    2023-10-31 [1] Bioconductor
#  lobstr                 * 1.1.2     2022-06-22 [1] CRAN (R 4.3.0)
#  locfit                   1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
#  magick                   2.8.1     2023-10-22 [1] CRAN (R 4.3.2)
#  magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  maps                     3.4.1.1   2023-11-03 [1] CRAN (R 4.3.2)
#  Matrix                 * 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
#  MatrixGenerics         * 1.14.0    2023-10-24 [1] Bioconductor
#  matrixStats            * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
#  memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
#  metapod                  1.10.0    2023-10-24 [1] Bioconductor
#  mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
#  munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
#  paletteer                1.6.0     2024-01-21 [1] CRAN (R 4.3.2)
#  pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  plotly                   4.10.3    2023-10-21 [1] CRAN (R 4.3.2)
#  png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.2)
#  prettyunits              1.2.0     2023-09-24 [1] CRAN (R 4.3.1)
#  promises                 1.2.1     2023-08-10 [1] CRAN (R 4.3.1)
#  purrr                    1.0.2     2023-08-10 [1] CRAN (R 4.3.1)
#  R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
#  rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                     1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
#  RCurl                    1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
#  rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.3.0)
#  restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.2)
#  rhdf5                  * 2.46.1    2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
#  rhdf5filters             1.14.1    2023-11-06 [1] Bioconductor
#  Rhdf5lib                 1.24.1    2023-12-11 [1] Bioconductor 3.18 (R 4.3.2)
#  rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.2)
#  rlang                    1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
#  rprojroot                2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
#  Rsamtools                2.18.0    2023-10-24 [1] Bioconductor
#  RSQLite                  2.3.4     2023-12-08 [1] CRAN (R 4.3.2)
#  rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.3.2)
#  rtracklayer              1.62.0    2023-10-24 [1] Bioconductor
#  S4Arrays               * 1.2.0     2023-10-24 [1] Bioconductor
#  S4Vectors              * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
#  sass                     0.4.8     2023-12-06 [1] CRAN (R 4.3.2)
#  ScaledMatrix             1.10.0    2023-10-24 [1] Bioconductor
#  scales                   1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
#  scater                   1.30.1    2023-11-16 [1] Bioconductor
#  scran                    1.30.0    2023-10-24 [1] Bioconductor
#  scuttle                  1.12.0    2023-10-24 [1] Bioconductor
#  sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
#  shiny                    1.8.0     2023-11-17 [1] CRAN (R 4.3.2)
#  shinyWidgets             0.8.0     2023-08-30 [1] CRAN (R 4.3.2)
#  SingleCellExperiment   * 1.24.0    2023-10-24 [1] Bioconductor
#  spam                     2.10-0    2023-10-23 [1] CRAN (R 4.3.2)
#  SparseArray            * 1.2.2     2023-11-07 [1] Bioconductor
#  sparseMatrixStats        1.14.0    2023-10-24 [1] Bioconductor
#  SpatialExperiment      * 1.12.0    2023-10-24 [1] Bioconductor
#  spatialLIBD            * 1.17.8    2024-09-11 [1] Github (LieberInstitute/spatialLIBD@c82c789)
#  statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
#  SummarizedExperiment   * 1.32.0    2023-10-24 [1] Bioconductor
#  tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
#  tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
#  utf8                     1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
#  vctrs                    0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
#  vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.3.2)
#  viridis                  0.6.4     2023-07-22 [1] CRAN (R 4.3.2)
#  viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
#  withr                    2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
#  XML                      3.99-0.16 2023-11-29 [1] CRAN (R 4.3.2)
#  xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                  0.42.0    2023-10-24 [1] Bioconductor
#  yaml                     2.3.8     2023-12-11 [1] CRAN (R 4.3.2)
#  zlibbioc                 1.48.0    2023-10-24 [1] Bioconductor

#  [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
