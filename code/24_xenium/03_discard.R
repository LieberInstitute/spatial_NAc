library(SpatialExperiment)
library(SpotSweeper)
library(sessioninfo)
library(ggplot2)
library(scater)
library(escheR)
library(here)

#Read in the post QC object 
spe <- readRDS(here("processed-data","xenium","SPEs","spe_splitQC_Global_prediscard.Rds"))

spe


#Remove cells based on cell area sum and detected
spe$discard <- spe$detected_outlier | spe$sum_outlier | spe$cell_area_outlier

table(spe$discard)

#violin plots of detected genes, sum UMIs. and cell area
#sum
sum_discard <- plotColData(object = spe,y = "sum",x = "Sample",colour_by = "discard") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = sum_discard,
       filename = here("plots","xenium","01_buildSPE","QC","sum_discard_violin_plot.png"))

#detected
detected_discard <- plotColData(object = spe,y = "detected",x = "Sample",colour_by = "discard") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = detected_discard,
       filename = here("plots","xenium","01_buildSPE","QC","detected_discard_violin_plot.png"))

#cell area
cell_area_discard <- plotColData(object = spe,y = "cell_area",x = "Sample",colour_by = "discard") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = cell_area_discard,
       filename = here("plots","xenium","01_buildSPE","QC","cell_area_discard_violin_plot.png"))

#Plot discard on tissue
for(sample in unique(spe$Sample)){
  print(sample)
  sub_spe <- spe[,spe$Sample == sample]
  p <- make_escheR(sub_spe) |> 
    add_fill("detected",point_size = 2) |>
    add_ground("discard",point_size = 2,stroke = 1) +
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "transparent"))
  ggsave(plot = p,
         filename = here("plots","xenium","01_buildSPE","QC",paste0(sample,"_Discard.png")),
         height = 18, width = 18)
}

#Remove the discarded cells
spe <- spe[,!spe$discard]

spe


#Save spe with QC information 
message(paste0("Saving clean spe - ",Sys.time()))
saveRDS(spe,here("processed-data","xenium","SPEs","spe_clean.Rds"))


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
