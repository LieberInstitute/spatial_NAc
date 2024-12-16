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

opt <- list()
opt$gene_selection_strategy <- "all_genes"
opt$data <- "human_NAc"
res_dir <- here("processed-data", "16_transfer_learning", "02_target_projections", opt$data)
plot_dir <- here("plots", "16_transfer_learning", "02_target_projections", opt$data)

spe <- readRDS(file = file.path(res_dir, paste0("spe_NMF_", opt$gene_selection_strategy, ".rds")))

# Specify the number of non-zero spots
if(opt$data == "human_NAc"){
    nFactors <- 83
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
pdf(file.path(plot_dir, paste0("num_nonzero_spots_", opt$gene_selection_strategy, ".pdf")), height = 6, width = 4)
ggplot(res_df[res_df$NMF %in% paste0("nmf", c(1:20)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
ggplot(res_df[res_df$NMF %in% paste0("nmf", c(21:39)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
ggplot(res_df[res_df$NMF %in% paste0("nmf", c(41:55)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
ggplot(res_df[res_df$NMF %in% paste0("nmf", c(56:83)), ], aes(x = NMF, y = num_spots, fill = NMF)) + coord_flip() + geom_bar(stat = "identity") + geom_hline(yintercept = 1000) + theme_classic() + theme(legend.position = "none")
dev.off()

res_df_final <- res_df[res_df$num_spots > 200, ]
dim(res_df_final)

pdf(file.path(plot_dir, paste0("distribution_nSpots_nonzero_coeff_", opt$gene_selection_strategy, ".pdf")), height = 3, width = 6)
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
plot_list <- lapply(res_df_final$NMF, function(n){
    plot_donor_list <- lapply(sample_order, function(isample){
        cat(isample, "\n")
         vis_gene(spe, sampleid = isample, geneid = n, is_stitched = TRUE) + ggtitle(isample)
    })
    plot_donor_list
})
names(plot_list) <- res_df_final$NMF

for(i in c(1:length(plot_list))){
    pdf(file.path(plot_dir, paste0(names(plot_list)[i], ".pdf")))
    print(plot_list[[i]])
    dev.off()
}

sessionInfo()
