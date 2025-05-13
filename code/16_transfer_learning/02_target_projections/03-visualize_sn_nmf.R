####snRNA-seq NMF pattern projection to Visium
library(tidyverse)
library(RcppML)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(here)
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


theme_Publication <- function(base_size=14, base_family="sans") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.box = "vetical",
               legend.key.size= unit(0.5, "cm"),
               #legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
       ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
      
}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#f87f01","#7fc97f","#ef3b2c","#feca01","#a6cee3","#fb9a99","#984ea3","#8C591D")), ...)
      
}

spec <- matrix(
    c("data", "d", 1, "character", "Specify the input dataset"
    ),
    byrow = TRUE, ncol = 5
)
opt <- getopt(spec)

#opt <- list()
#opt$data <- "human_NAc"
print(opt$data)
res_dir <- here("processed-data", "16_transfer_learning", "02_target_projections", opt$data)
plot_dir <- here("plots", "16_transfer_learning", "02_target_projections", opt$data)

spe <- readRDS(file = file.path(res_dir, paste0("spe_NMF.rds")))

# Specify the number of non-zero spots
if(opt$data == "human_NAc"){
    nFactors <- 66
}else{
    if(opt$data == "rat_case_control_acute"){
        nFactors <- 43
    }else{
        if(opt$data == "rat_case_control_repeated"){
            nFactors <- 40
        }else{
            if(opt$data == "rat_case_control_morphine_acute"){
                nFactors <- 46
            }else{
                if(opt$data == "rat_case_control_morphine_repeated"){
                    nFactors <- 36
                }
            }
        }
    }
}

nmf_factors <- paste0("nmf", c(1:nFactors))
nSpot_nonzero_coeff <- c()
for(i in c(1:length(nmf_factors))){
    nSpot_nonzero_coeff[i] <- sum(spe[[nmf_factors[i]]] > 0)
}

nSpot_nonzero_coeff_by_sample <- lapply(c(1:nFactors), function(i){
    table(spe$sample_id[spe[[nmf_factors[i]]] > 0])
})
nSpot_nonzero_coeff_by_sample <- do.call(rbind, nSpot_nonzero_coeff_by_sample)
nSpot_nonzero_coeff_by_sample <- data.frame(nSpot_nonzero_coeff_by_sample)
nSpot_nonzero_coeff_by_sample$NMF <- c(1:nFactors)
nSpot_nonzero_coeff_by_sample <- reshape2::melt(nSpot_nonzero_coeff_by_sample, id.vars = "NMF")

res_df <- data.frame("NMF" = nmf_factors, "num_spots" = nSpot_nonzero_coeff)
res_df$NMF <- factor(res_df$NMF, levels = nmf_factors)

if(opt$data == "human_NAc"){
    pdf(file.path(plot_dir, paste0("num_nonzero_spots.pdf")), height = 6, width = 4)
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:39)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(40:66)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    dev.off()
}
if(opt$data == "rat_case_control_acute"){
    pdf(file.path(plot_dir, paste0("num_nonzero_spots.pdf")), height = 6, width = 4)
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:43)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    dev.off()
} 
if(opt$data == "rat_case_control_repeated"){
    pdf(file.path(plot_dir, paste0("num_nonzero_spots.pdf")), height = 6, width = 4)
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:40)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    dev.off()
}
if(opt$data == "rat_case_control_morphine_acute"){
    pdf(file.path(plot_dir, paste0("num_nonzero_spots.pdf")), height = 6, width = 4)
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:46)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    dev.off()
}
if(opt$data == "rat_case_control_morphine_repeated"){
    pdf(file.path(plot_dir, paste0("num_nonzero_spots.pdf")), height = 6, width = 4)
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:36)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
    dev.off()
}

res_df_final <- res_df[res_df$num_spots > 200, ]

pdf(file.path(plot_dir, paste0("distribution_nSpots_nonzero_coeff.pdf")), height = 3, width = 6)
ggplot(res_df_final, aes(x = num_spots)) + geom_histogram() + theme_classic() + xlab("Number of spots w/ non-zero coefficient") + ylab("Number of factors")
dev.off()


# Visualization of patterns
spe_Br2743 <- mirrorObject(spe, sample_id = "Br2743", image_id = "lowres", axis = "v")
spe_Br8492 <- mirrorObject(spe, sample_id = "Br8492", image_id = "lowres", axis = "v")
spe_Br8325 <- mirrorObject(spe, sample_id = "Br8325", image_id = "lowres", axis = "v")
spe_Br3942 <- mirrorObject(spe, sample_id = "Br3942", image_id = "lowres", axis = "v")

spe <- spe[ ,!spe$sample_id %in% c("Br2743", "Br8492", "Br8325", "Br3942")]
spe <- cbind(spe, spe_Br2743)
spe <- cbind(spe, spe_Br8492)
spe <- cbind(spe, spe_Br8325)
spe <- cbind(spe, spe_Br3942)
sample_order <- c("Br2743", "Br6432", "Br6423", "Br2720", "Br6471", "Br6522","Br8492", "Br8325", "Br8667", "Br3942")

# Add clustering data
safe_colorblind_palette <- c("#E7298A", "#A6761D","#D95F02" , "#E6AB02",  "#66A61E","#1B9E77", "#7570B3","#666666")
clustDir <- here("processed-data", "07_spatial_domains", "01_precast", "nnSVG_precast", "final_clusters")
clusters_df <- read.csv(file.path(clustDir, "precast_clusters.csv"))
spe <- spe[ ,colnames(spe) %in% clusters_df$key]
clusters_df <- clusters_df[match(colnames(spe), clusters_df$key), ]
spe$PRECAST_clusters <- clusters_df$cluster
spe$PRECAST_clusters <- factor(spe$PRECAST_clusters, 
levels = c("D1 islands", "Endothelial/Ependymal", "Excitatory", "Inhibitory", "MSN 1", "MSN 2", "MSN 3", "WM"))


for(n in levels(res_df$NMF)){
     plot_donor_list <- lapply(sample_order, function(isample){
        cat(isample, "\n")
         vis_gene(spe, sampleid = isample, geneid = n, is_stitched = TRUE) + ggtitle(isample)
    })
    pdf(file.path(plot_dir, paste0(n, ".pdf")))
    print(plot_donor_list)
    dev.off()
}

for(n in levels(res_df$NMF)){
    plot_donor_list <- lapply(sample_order, function(isample){
        cat(isample, "\n")
        spe_single <- spe[ ,spe$sample_id == isample]
        spe_single <- spe_single[ ,!spe_single$exclude_overlapping]
        p <-  make_escheR(spe_single) |>
        add_fill(var = as.character(n)) |>
        add_ground(var = "PRECAST_clusters", stroke = 0.5) + scale_color_manual(values = safe_colorblind_palette) +
        scale_fill_gradient(low = "white", high = "black")+
        guides(color=guide_legend(title="PRECAST clusters")) +theme(plot.title = element_text(size = 35, face = "bold"),
            legend.title=element_text(size=20), 
            legend.text=element_text(size=18)) + ggtitle(isample)
        p
    })
    pdf(file.path(plot_dir, paste0(n, "_escheR.pdf")), height = 11, width = 18)
    print(plot_donor_list)
    dev.off()
}


barplot_list <- list()
for(i in res_df_final$NMF){
    spe_subset <- spe[ ,spe[[i]] > 0]
    df <- colData(spe_subset)
    df <- data.frame(df)
    barplot_list[[as.character(i)]] <- ggplot(df, aes_string(x = "PRECAST_clusters", y = as.character(i), fill = "PRECAST_clusters")) + geom_violin() + xlab("Spatial domains") + ylab(toupper(as.character(i))) + scale_fill_manual(values = safe_colorblind_palette) + theme_Publication() + theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + guides(fill=guide_legend(title="Spatial domains")) + ggtitle(toupper(as.character(i)))
}

pdf(file.path(plot_dir, "NMF_coeff_by_spatial_domains.pdf"), width = 10, height = 5)
for(i in c(1:length(barplot_list))){
    print(barplot_list[[i]])
}
dev.off()

sessionInfo()
