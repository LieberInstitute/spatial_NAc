#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#code modified from: https://github.com/LieberInstitute/spatialdACC/blob/main/code/snRNA-seq/01_QC/01_empty_drops.R

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(here)
library(sessioninfo)
library(ggplot2)
library(tidyverse)

## get sample i
sample_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#### Load & Subset raw data ####
sce <- readRDS(here("processed-data","12_snRNA","sce_raw.Rds"))

#Print some info about sce object
as.data.frame(unique(colData(sce)[,c("Sample","Raw_data_path")]))

#Make sure rownames of the colData match the column names on the count matrix. Stop, if not. 
identical(rownames(colData(sce)),colnames(sce))
stopifnot(identical(rownames(colData(sce)),colnames(sce)))

samples <- unique(sce$Sample)
sample_run <- samples[[sample_i]]
message("Running Sample: ", sample_run, " (", sample_i, "/", length(samples), ")")

sce <- sce[, sce$Sample == sample_run]
message("ncol:", ncol(sce))

#### Run barcodeRanks to find knee ####

bcRanks <- barcodeRanks(sce, fit.bounds = c(10, 1e3))

knee_highest <- metadata(bcRanks)$knee - 200
message(
  "'First knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_highest =", knee_highest
)

knee_higher <- metadata(bcRanks)$knee - 100
message(
  "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_higher =", knee_higher
)

message(
  "'Third knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee =", metadata(bcRanks)$knee
)

knee_lower <- metadata(bcRanks)$knee + 100
message(
  "'Fourth knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_lower =", knee_lower
)

knee_lowest <- metadata(bcRanks)$knee + 200
message(
  "'Fifth knee point' = ", metadata(bcRanks)$knee, "\n",
  "knee_lowest =", knee_lowest
)

#### Run emptyDrops w/ knee + 100 ####
set.seed(1234)
message("Starting emptyDrops")
Sys.time()
e.out <- DropletUtils::emptyDrops(
  sce,
  niters = 30000,
  lower = knee_lower
)
message("Done - saving data")
Sys.time()

save(e.out,
     file = here("processed-data","12_snRNA","droplet_scores",paste0(sample_run,"_droplet_scores.Rdata")))

#### QC Plots ####
message("QC check")
FDR_cutoff <- 0.001
addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))

n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
message(n_cell_anno)

my_theme <- theme_bw() +
  theme(text = element_text(size = 15))

droplet_elbow_plot <- as.data.frame(bcRanks) %>%
  add_column(FDR = e.out$FDR) %>%
  ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
  annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Knee", vjust = -1, color = "gray") +
  geom_hline(yintercept = knee_highest, linetype = "dashed") +
  annotate("text", x = 10, y = knee_highest, label = "Knee est 'highest'") +
  geom_hline(yintercept = knee_higher, linetype = "dashed") +
  annotate("text", x = 10, y = knee_higher, label = "Knee est 'higher'") +
  geom_hline(yintercept = knee_lower, linetype = "dashed") +
  annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'") +
  geom_hline(yintercept = knee_lowest, linetype = "dashed") +
  annotate("text", x = 10, y = knee_lowest, label = "Knee est 'lowest'") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "Barcode Rank",
    y = "Total UMIs",
    title = paste("Sample:", sample_run),
    subtitle = n_cell_anno,
    color = paste("FDR <", FDR_cutoff)
  ) +
  my_theme +
  theme(legend.position = "bottom")

ggsave(droplet_elbow_plot, 
       filename = here("plots","12_snRNA","droplet_scores",paste0(sample_run,"_droplet_qc.png")))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
