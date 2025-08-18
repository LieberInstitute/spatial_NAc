# Load libraries
library(here)
library(ggplot2)
library(tidyverse)

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.85)
    as.numeric(column >= top_10_percent_threshold)
}

res_dir <- here("processed-data", "19_sLDSC", "10x_visium", "input_files")
plot_dir <- here("plots", "19_sLDSC", "10x_visium")
dat <- read.table(file.path(res_dir, "10x_visium_aggregated_cpm.tsv"),header=T)
dat.norm <- t(apply(dat, 1, function(x){x/sum(x)}))
res <- apply(dat.norm, 2, is_top_10_percent)
rownames(res) <- rownames(dat.norm)
write.csv(res,file.path(res_dir, "10x_visium_score.csv"))

# Check for the overlap between known disease genes
openTargets_Dir <- here("processed-data", "OpenTargets")
depression_genes <- data.frame(read_tsv(file.path(openTargets_Dir, "depression_open_targets.tsv")))
substance_abuse_genes <- data.frame(read_tsv(file.path(openTargets_Dir, "substance_abuse_open_targets.tsv")))
opioid_dependence_genes <- data.frame(read_tsv(file.path(openTargets_Dir, "opioid_dependence_open_targets.tsv")))
schizophrenia_genes <- data.frame(read_tsv(file.path(openTargets_Dir, "schizophrenia_open_targets.tsv")))
bipolar_genes <- data.frame(read_tsv(file.path(openTargets_Dir, "bipolar_open_targets.tsv")))

depression_genes <- depression_genes[depression_genes$globalScore > 0.1, ]
substance_abuse_genes <- substance_abuse_genes[substance_abuse_genes$globalScore > 0.1, ]
opioid_dependence_genes <- opioid_dependence_genes[opioid_dependence_genes$globalScore > 0.1, ]
schizophrenia_genes <- schizophrenia_genes[schizophrenia_genes$globalScore > 0.1, ]
bipolar_genes <- bipolar_genes[bipolar_genes$globalScore > 0.1, ]

get_intersection <- function(df1, df2, icell){
    genes_select <- rownames(df1)[df1[ ,icell] == 1]
    length(intersect(genes_select, df2$symbol))
}
depression_intersect <- c()
substance_abuse_intersect <- c()
opioid_dependence_intersect <- c()
bipolar_intersect <- c()
schizophrenia_intersect <- c()
for(i in colnames(res)){
   depression_intersect[i] <- get_intersection(res, depression_genes, i)
   substance_abuse_intersect[i] <- get_intersection(res, substance_abuse_genes, i)
   opioid_dependence_intersect[i] <- get_intersection(res, opioid_dependence_genes, i)
   bipolar_intersect[i] <- get_intersection(res, bipolar_genes, i)
   schizophrenia_intersect[i] <- get_intersection(res, schizophrenia_genes, i)
}

depression_intersect <- data.frame(depression_intersect)
substance_abuse_intersect <- data.frame(substance_abuse_intersect)
opioid_dependence_intersect <- data.frame(opioid_dependence_intersect)
bipolar_intersect <- data.frame(bipolar_intersect)
schizophrenia_intersect <- data.frame(schizophrenia_intersect)

depression_intersect$cellTypes <- rownames(depression_intersect)
substance_abuse_intersect$cellTypes <- rownames(substance_abuse_intersect)
bipolar_intersect$cellTypes <- rownames(bipolar_intersect)
opioid_dependence_intersect$cellTypes <- rownames(opioid_dependence_intersect)
schizophrenia_intersect$cellTypes <- rownames(schizophrenia_intersect)

colnames(depression_intersect) <- c("nGenes", "cellType")
colnames(substance_abuse_intersect) <- c("nGenes", "cellType")
colnames(bipolar_intersect) <- c("nGenes", "cellType")
colnames(opioid_dependence_intersect) <- c("nGenes", "cellType")
colnames(schizophrenia_intersect) <- c("nGenes", "cellType")

plot_list <- list()
plot_list[[1]] <- ggplot(depression_intersect, aes(x = cellType, y = nGenes)) + geom_bar(stat = "identity") + ggtitle("Depression") + xlab("") + ylab("nGenes") + coord_flip()
plot_list[[2]] <- ggplot(substance_abuse_intersect, aes(x = cellType, y = nGenes)) + geom_bar(stat = "identity") + ggtitle("Substance abuse") + xlab("") + ylab("nGenes") + coord_flip()
plot_list[[3]] <- ggplot(opioid_dependence_intersect, aes(x = cellType, y = nGenes)) + geom_bar(stat = "identity") + ggtitle("Opioid dependence") + xlab("") + ylab("nGenes") + coord_flip()
plot_list[[4]] <- ggplot(bipolar_intersect, aes(x = cellType, y = nGenes)) + geom_bar(stat = "identity") + ggtitle("Bipolar") + xlab("") + ylab("nGenes") + coord_flip()
plot_list[[5]] <- ggplot(schizophrenia_intersect, aes(x = cellType, y = nGenes)) + geom_bar(stat = "identity") + ggtitle("Schizophrenia") + xlab("") + ylab("nGenes") + coord_flip()

pdf(file.path(plot_dir, "disease_genes_overlap.pdf"), width = 4, height = 7)
for(i in c(1:length(plot_list))){
    print(plot_list[[i]])
}
dev.off()

overlap_matrix <- matrix(NA, nrow = dim(res)[2], ncol = dim(res)[2])
for(i in c(1:dim(res)[2])){
    for(j in c(1:dim(res)[2])){
        overlap_matrix[i,j] <- length(intersect(rownames(res)[res[ ,i] == 1], rownames(res)[res[ ,j] == 1]))
    }
}
overlap_matrix <- t(apply(overlap_matrix, 1, function(x){x/max(x)}))
rownames(overlap_matrix) <- colnames(res)
colnames(overlap_matrix) <- colnames(res)

library(corrplot)
rownames(overlap_matrix) <- gsub("_", " ", rownames(overlap_matrix))
colnames(overlap_matrix) <- gsub("_", " ", colnames(overlap_matrix))
pdf(file.path(plot_dir, "overlap_between_categories.pdf"), width = 6, height = 5)
corrplot(overlap_matrix, type = "lower", order = "hclust", method = "shade", is.corr = FALSE, tl.cex = 0.8, tl.col = "black")
dev.off()






