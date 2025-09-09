#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#module load r_nac

library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(ggplot2)
library(scater)
library(dplyr)
library(scran)
library(here)

sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))


#PPP1R1B featureplot
PPP1R1B_Feature <- plotReducedDim(object = sce,
                               dimred = "tSNE_HARMONY",
                               color_by = "PPP1R1B",
                               swap_rownames = "gene_name") +
  scale_color_gradientn(colours = c("lightgrey","red")) +
  theme(axis.line = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

#with legend
ggsave(plot = PPP1R1B_Feature,
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","PPP1R1B_FeaturePlot_legend.pdf"),
       height = 8,
       width = 8)
#without legend
ggsave(plot = PPP1R1B_Feature + theme(legend.position = "none"),
       filename = here("plots","12_snRNA","Supplementary","Additional_Feature_Plots","PPP1R1B_FeaturePlot_nolegend.png"),
       height = 8,
       width = 8)

sessionInfo()
