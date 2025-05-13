####code for correlating categorical technical variables with NMF patterns
library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(here)

opt <- list()
opt$data <- "human_NAc"
print(opt$data)

res_dir <- here::here("processed-data", "16_transfer_learning","01_process_reference", "RCppML", opt$data)
plot_dir <- here::here("plots", "16_transfer_learning","01_process_reference", "RCppML", opt$data)

x <- readRDS(file = file.path(res_dir,paste0("nmf_results.rds")))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# Prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h

#custom plotting function
gseaPlot <- function(nmf, nmf_terms) {
  file1 = grep(nmf, list.files(res_dir), value=T)
  results = read.csv(paste0(res_dir,"/", file1), row.names=1)
  tmp = arrange(results[results$pathway %in% nmf_terms,], desc(padj))
  #pull rownames of nonzero
  non0.nmf = rownames(loads)[loads[,nmf]>0]
  #convert to entrez id to match with reactome entires
  non0.id = mapIds(org.Hs.eg.db, keys=non0.nmf, keytype="ENSEMBL", column="ENTREZID", multiVals = "first")
  #match to gene names
  names(non0.id) = non0.nmf
  #remove gene names without ID
  non0.id = non0.id[!is.na(non0.id)]
  #pull nmf values of mapped genes
  nmf.stats = loads[names(non0.id),nmf]
  #convert names to entrez ID
  names(nmf.stats) = non0.id
  #place in rank order
  nmf.stats = sort(nmf.stats, decreasing=T)
  #make into df
  nmf.df = do.call(rbind, lapply(1:nrow(tmp), function(x) {
    term1 = tmp$pathway[x]
    as.data.frame(list("id"=names(nmf.stats),
                       "weight"=nmf.stats,
                       "rank"=1:length(nmf.stats),
                       "term"=term1,
                       "present"=names(nmf.stats) %in% x.h[[term1]],
                       y=x-1, yend=x))
  }))
  y_labels= sapply(1:nrow(tmp), function(x) paste0("NES= ",round(tmp[x,"NES"],2),
                                                   "\n(",tmp[x,"size"],"/ ",length(x.h[[tmp[x,"pathway"]]]),")\n",
                                                   "padj= ",format(tmp[x,"padj"], scientific=T, digits=2)))
  ggplot(filter(nmf.df, present==T), aes(x=rank, xend=rank, y=y, yend=yend, color=term))+
    geom_segment()+scale_x_continuous(limits=c(0,max(nmf.df$rank)), expand=c(0,0))+
    scale_y_continuous(breaks=c(0.5, 1.5, 2.5, 3.5, 4.5), labels=y_labels, expand=c(0,0))+
    theme_minimal()+ggtitle(nmf)+
    theme(axis.title.y=element_blank(), legend.title=element_blank(),
          legend.position="bottom", legend.direction = "vertical",
          panel.grid.minor=element_blank(), panel.grid.major=element_blank())
}

nmf10.terms = c("Signaling by GPCR","G alpha (i) signalling events","G alpha (q) signalling events", "Opioid Signalling", "Class A/1 (Rhodopsin-like receptors)")
p1 <- gseaPlot("nmf10", nmf10.terms)

nmf7.terms = c("Signaling by GPCR", "Transmission across Chemical Synapses", "Opioid Signalling", "Activation of NMDA receptors and postsynaptic events", "Neurexins and neuroligins")
p2 <- gseaPlot("nmf7", nmf7.terms)

nmf4.terms = c("Signaling by GPCR", "Transmission across Chemical Synapses", "Protein-protein interactions at synapses", "G alpha (i) signalling events", "Hedgehog ligand biogenesis")
p3 <- gseaPlot("nmf4", nmf4.terms)

nmf3.terms = c("Signaling by GPCR", "Transmission across Chemical Synapses", "Opioid Signalling", "Sensory Perception", "Neurexins and neuroligins")
p4 <- gseaPlot("nmf3", nmf3.terms)

ggsave(file.path(plot_dir, "Figure3_custom-gsea-plots.pdf"), gridExtra::grid.arrange(p4, p3, p2, p1, ncol=1),
       bg="white", height=20, width=12, units="in")