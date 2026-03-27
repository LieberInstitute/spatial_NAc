#Goal: Calculate size factors based on nucleus and cell area
library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)

spe <- readRDS(here("processed-data","xenium","SPEs","spe_clean.Rds"))
spe

#Based on Atta et al. 2024, generate non-count based scaling factors
#Calculate both cell area and nucleus area values
spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area)

nuc_area_hist <- ggplot(colData(spe),aes(x = nucleus_area.sf)) +
  geom_histogram(bins = 50) +
  labs(x = "Scaling Factor",
       y = "Frequency",
       title = "Nucleus area") +
  theme_bw() +
  theme(plot.title = element_text(hjust= 0.5))
ggsave(plot = nuc_area_hist,filename = here("plots","xenium","nucleus_area_scaling_hist.pdf"))

#Lines 58-61 straight from: https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/Xenium/03_quality_control/01_perCellQC.R
# normalize the counts by the nucleus and cell area scaling factors
assay(spe, "nucleus_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$nucleus_area.sf, transform="log", assay.type="counts")


#Plot genes to help identify WM tracts and NAc proper
spe$MOBP <- assay(spe,"nucleus_normcounts")["MOBP",]
spe$ST18 <- assay(spe,"nucleus_normcounts")["ST18",]
spe$DRD1 <- assay(spe,"nucleus_normcounts")["DRD1",]
spe$RXFP1 <- assay(spe,"nucleus_normcounts")["RXFP1",]
spe$OPRM1 <- assay(spe,"nucleus_normcounts")["OPRM1",]
spe$PROK2 <- assay(spe,"nucleus_normcounts")["PROK2",]
spe$GABRQ <- assay(spe,"nucleus_normcounts")["GABRQ",]

for(sample in unique(spe$Sample)){
  print(sample)
  sub_spe <- spe[,spe$Sample == sample]
  for(gene in c("RXFP1","OPRM1","PROK2","GABRQ","MOBP","ST18","DRD1")){
    print(gene)
    p1 <- make_escheR(sub_spe) |>
      add_fill(gene) +
      scale_fill_gradientn(colors = c("grey","red"))
    ggsave(plot = p1,
           filename = here("plots","xenium","Expression",
                           paste0(sample,"_",gene,".png")),
           height = 18,width = 18)
  }
}

#Save spe with normalized counts
saveRDS(spe,here("processed-data","xenium","SPEs","spe_NormCounts.Rds"))


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
