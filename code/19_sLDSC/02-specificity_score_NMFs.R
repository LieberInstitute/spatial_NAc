library(here)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.9)
    as.numeric(column >= top_10_percent_threshold)
}

res_dir <- here("processed-data", "19_sLDSC", "human_NAc_NMF", "input_files")
plot_dir <- here("plots", "19_sLDSC", "human_NAc_NMF")
dat <- read.table(file.path(res_dir, "NMF_correlations.tsv"),header=T)
res <- apply(dat, 2, is_top_10_percent)
rownames(res) <- rownames(dat)
write.csv(res,file.path(res_dir, "NMFs_score.csv"))

plot_df1 <- dat
plot_df1$gene <- rownames(plot_df1)
plot_df2 <- data.frame(res) 
plot_df2$gene <- rownames(plot_df2)
plot_df1 <- reshape2::melt(plot_df1, id.vars = "gene")
plot_df2 <- reshape2::melt(plot_df2, id.vars = "gene")
colnames(plot_df1) <- c("gene", "nmf", "gene_corr")
colnames(plot_df2) <- c("gene", "nmf", "select")
plot_df <- merge(plot_df1, plot_df2, by=c("gene","nmf"))

plot_df$select <- factor(plot_df$select, levels = c(0, 1))
pdf(file.path(plot_dir, "correlation_distribution_by_NMF.pdf"), height = 12, width = 5)
ggplot(plot_df, aes(x = nmf, y = gene_corr, fill = select)) + geom_boxplot(outlier.size = 0.5) + theme_classic() + coord_flip() + xlab("") + ylab("Gene correlation") +
geom_hline(yintercept = 0, color = "red", linetype="dashed")
dev.off()


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

intersect_df <- data.frame("Depression" = depression_intersect, 
                           "Substance_abuse" = substance_abuse_intersect,
                           "Opioid_dependence" = opioid_dependence_intersect, 
                           "Bipolar" = bipolar_intersect, 
                           "Schizophrenia" = schizophrenia_intersect)
intersect_df$NMF <- rownames(intersect_df)
intersect_df <- reshape2::melt(intersect_df, id.vars = "NMF")
intersect_df$NMF <- factor(intersect_df$NMF, levels = paste0("nmf", c(1:66)))

pdf(file.path(plot_dir, "disease_genes_overlap.pdf"), width = 8, height = 10)
ggplot(intersect_df, aes(x = variable, y = NMF, fill = value)) + geom_tile() +
  scale_fill_gradient(low="white", high="blue") + xlab("") + ylab("") + theme_minimal()
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
pdf(file.path(plot_dir, "overlap_between_categories.pdf"), width = 10, height = 10)
corrplot(overlap_matrix, type = "lower", order = "hclust", method = "shade", is.corr = FALSE, tl.cex = 0.8, tl.col = "black")
dev.off()



