rm(list = ls())
setwd('/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/')
library(biomaRt)
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)
library(here)
library(tidyr)
library(forcats)
library(ggplot2)
library(dplyr)
library(patchwork)


##load the data
mch<-readH5AD(file=here::here('processed-data','25_revisions',
                              '03-projections_retroSeq','rs2_mch_matrix.h5ad'))

##get gene names
mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", values = rownames(mch),
              mart = mart)
symbs <- symb$mgi_symbol[match(rownames(mch), symb$ensembl_gene_id, nomatch = NA)]
sum(is.na(symbs))
#[1] 528

rowData(mch)$gene_name<-symbs
rowData(mch)$start<-NULL
rowData(mch)$end<-NULL

##drop genes with no names (need these for ortholog matching)
mch<-mch[!is.na(rowData(mch)$gene_name),]
rownames(mch)<-rowData(mch)$gene_name

# Translate from one mchcies to the other using the orthology
orthology<-read.csv(file=here::here('processed-data','25_revisions','03-projections_retroSeq',
                                    'human_mouse_orthologs.csv'))
names <- orthology[orthology$Column3 %in% rownames(mch),]

names <- names[match(rownames(mch), names$Column3),]

setdiff(names$Column3, rownames(mch))

sum(is.na(names$Column1))
# [1] 12420

rownames(mch) <- names$Column1
dim(mch)
# [1] 31680  262

# remove rownames that are NA
mch <- mch[!is.na(rownames(mch)),]
dim(mch)
# [1] 19260  262

# there are some repeated rownames in mch
duplicated_rows <- duplicated(rownames(mch))
mch <- mch[!duplicated_rows, ]
dim(mch)
# [1] 17458  262

##load nmf patterns
nmf_dir <- here::here("processed-data", "16_transfer_learning", "01_process_reference", "RCppML", "human_NAc")
sce_dir <- here::here("processed-data", "12_snRNA")
sce <- readRDS(file = file.path(sce_dir, "sce_CellType_noresiduals.Rds"))
geneData <- rowData(sce)

nmf_results <- readRDS(file.path(nmf_dir, "nmf_results.rds"))
loadings <- nmf_results@w

no_expr <- which(rowSums(loadings) == 0)
loadings <- loadings[-no_expr, ]
dim(loadings)
# 14002    66

geneData <- geneData[geneData$gene_id %in% rownames(loadings), ]
geneData <- geneData[match(rownames(loadings), geneData$gene_id), ]

# Convert loadings rownames to symbols
rownames(loadings) <- geneData$gene_name
set.seed(101)
i<-intersect(rownames(mch),rownames(loadings))
length(i)
# [1] 11628

loadings<-loadings[rownames(loadings) %in% i,]
mch<-mch[rownames(mch) %in% i,]


## projection
loadings<-loadings[match(rownames(mch),rownames(loadings)),]
proj<-project(loadings,assay(mch,'X'),L1=0)
proj<-t(proj)

proj<-apply(proj,2,function(x){x/sum(x)})

colData(mch)<-cbind(colData(mch),proj)

saveRDS(mch, file = here::here('processed-data', '25_revisions','03-projections_retroSeq', "mch_projection.rds"))

# Make a simple violin plot to compare
data<-as.data.frame(colData(mch))
# keep cols starting with "nmf"
nmf_cols <- grep("^nmf", colnames(data), value = TRUE)
data <- data[, c("Target",nmf_cols)]

data_clean <- data[, colSums(!is.na(data)) > 0]

# Run t-tests to check if any NMF is associated with Target
nmf_cols <- grep("^nmf", colnames(data_clean), value = TRUE)
results <- lapply(nmf_cols, function(col) {
  x <- data_clean[data_clean$Target == "HY", col]
  y <- data_clean[data_clean$Target == "VTA", col]
  
  # remove NA values
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  # only run test if both groups have variation
  if (length(unique(c(x, y))) > 1) {
    test <- t.test(x, y)
    
    data.frame(
      feature = col,
      mean_HY = mean(x),
      mean_VTA = mean(y),
      p_value = test$p.value
    )
  } else {
    NULL
  }
})

results_df <- do.call(rbind, results)
results_df$adj_p <- p.adjust(results_df$p_value, method = "BH")
sig_results <- results_df[results_df$adj_p < 0.05, ]

sig_features <- sig_results$feature

# optional: order NMFs by adjusted p-value
nmf_order <- sig_results %>%
  arrange(adj_p) %>%
  pull(feature)

#-----------------------------
# 1) compute proportion > 0
#-----------------------------
seed1 <- as.matrix(data_clean[, sig_features, drop = FALSE])
seed1 <- seed1 > 0

d1 <- cbind.data.frame(Target = data_clean$Target, seed1) %>%
  group_by(Target) %>%
  add_tally(name = "total") %>%
  group_by(Target, total) %>%
  summarise(across(all_of(sig_features), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(sig_features),
    names_to = "nmf",
    values_to = "n"
  ) %>%
  mutate(prop = n / total)

#-----------------------------
# 2) compute scaled average
#-----------------------------
seed2 <- scale(as.matrix(data_clean[, sig_features, drop = FALSE]))

d2 <- cbind.data.frame(Target = data_clean$Target, seed2) %>%
  group_by(Target) %>%
  summarise(across(all_of(sig_features), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(
    cols = all_of(sig_features),
    names_to = "nmf",
    values_to = "scaled.avg"
  )

#-----------------------------
# 3) merge
#-----------------------------
dot_df <- left_join(
  d1[, c("Target", "nmf", "prop")],
  d2[, c("Target", "nmf", "scaled.avg")],
  by = c("Target", "nmf")
)

dot_df$nmf_f <- factor(dot_df$nmf, levels = nmf_order)
dot_df$Target <- factor(dot_df$Target, levels = c("HY", "VTA"))

#-----------------------------
# 4) plot
#-----------------------------
plotDir <- here::here("plots", "25_revisions", "03-projections_retroSeq")
pdf(file.path(plotDir, "NMF_by_Target_dotplot.pdf"), width = 11, height = 4)

ggplot(dot_df, aes(x = nmf_f, y = Target, size = prop, color = scaled.avg)) +
  geom_point() +
  theme_bw() +
  scale_size(range = c(0.5, 6), name = "Proportion") +
  scale_color_viridis_c(option = "F", direction = -1, name = "Scaled Avg.") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  xlab("") +
  ylab("")

dev.off()

pdf(file.path(plotDir, "Target_counts_barplot.pdf"), width = 4, height = 4)

ggplot(data, aes(x = Target, fill = Target)) +
  geom_bar() + 
  theme_bw() +
  scale_fill_manual(values = c("HY" = "steelblue", "VTA" = "tomato")) +
  labs(y = "Count", x = "", title = "Cell Counts by Target") +
  theme(legend.position = "none")
dev.off()