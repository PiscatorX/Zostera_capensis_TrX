#!/usr/bin/env Rscript
rm(list=ls())
cat("\f")
library(EnhancedVolcano)
library(RColorBrewer)
library(BiocParallel)
library(tidyverse)
library(tximport)
library(pheatmap)
library(tximeta)
library(DESeq2)



register(MulticoreParam(4))

source("/home/drewx/Documents/Zostera_capensis_TrX/scripts/Xpression_funcs.R")

sampleTable <- read.table("/opt/DB_REF/Metadata/Z_capensis.Transcriptome2.tsv", sep = "\t", header = T)
quant_files <- file.path("/home/drewx/Documents/Zostera_capensis_TrX/scripts/data.enterprise/Salmon-Quant",paste0(sampleTable$Readname,"_rawlib.sf"))

names(quant_files) <- sampleTable$Readname

################################# Swissprot hits ###############################

(Evigene_hits = read.delim("data.enterprise/Blast_top_hits/EvigeneX.blastx.ids", sep="\t", header = F) %>%
                select("V1","V2") %>%
                `colnames<-`(c("contig","ACC"))) %>%
                 head()
dim(Evigene_hits)

(uniprot_ref <- read.delim("data/Blast_top_hits/uniprot_sprot.ref", sep="\t", header = F) %>% 
    `colnames<-`(c("ACC","name"))) %>%
     head()

(Uniprot_hits <- merge(Evigene_hits, uniprot_ref, all.x = T) %>%
                 select(contig, name)) %>%
                 head()

trX_data <- tximport(quant_files, type="salmon", tx2gene = Evigene_hits, txOut=F)

################################################################
#   trX_data <- tximport(quant_files, type="salmon", txOut=T)  #
################################################################

colnames(trX_data$counts)
################################# DESeq ########################################
rownames(sampleTable) <- colnames(trX_data$counts)
dds_dataset <- DESeqDataSetFromTximport(trX_data, sampleTable, ~Treatment)
dds_dexp <- DESeq(dds_dataset)
res <- results(dds_dexp, parallel = TRUE)
dds_vst <- vst(dds_dexp, blind=FALSE)
################################# PCA analysis #################################
mypca <- list()
mypca$pca_data <- plotPCA(dds_vst, intgroup ="Treatment", returnData=TRUE)
mypca$var_perc <- round(100 * attr(mypca$pca_data, "percentVar"))
mypca$pca_data$name <- factor(c("Control_rep1", "Control_rep2", "Control_rep3", "Heat_stressed_rep1", "Heat_stressed_rep2","Heat_stressed_rep3"))
px3 <- mypca_plot(mypca, label="name", plot_fname = "PCAplot.pdf")


########################################################### significant DE analysis ###########################################
#https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
  
signf_genes <- function(padj_cutoff){
  
    res_sig <- res%>%
               data.frame() %>%
               rownames_to_column(var="Transcript") %>%
               filter(padj < {padj_cutoff}) %>%
               filter(abs(log2FoldChange) > 1) %>%
               arrange(desc(log2FoldChange)) %>%
               as_tibble

return(res_sig)

}

(sig_genes1 <- signf_genes(padj_cutoff =  0.05))

res_dat <- res %>%
           data.frame %>%
           rownames_to_column(var="SwissProt")

write.table(res_dat, file = "./data.enterprise/Blast_top_hits/TrX_padj1.tsv", row.names = FALSE, quote = FALSE, sep ="\t", col.names = T)

(top_DEG_names <- sig_genes1 %>% pull(Transcript))

sig_genes1 <- sig_genes1 %>%
                rename(`Swiss-Prot`=Transcript)

#erge(


top_DEG_order <- order(rowMeans(counts(dds_dexp[top_DEG_names,], normalized=TRUE)), decreasing=TRUE)
  
########################################################### heatmap  ##################################################
  
colData(dds_vst)
  
metadata_col <- as.data.frame(colData(dds_vst)[,c("Rep_name")])
colnames(metadata_col) <-  c("Replicate")
metadata_row <- as.data.frame(colData(dds_vst)[,c("Rep_name")])
colnames(metadata_row) <-  c("Sample")
rownames(metadata_row) <- sampleTable$Readname
  
colnames(dds_vst) <- sampleTable$Rep_name
(data <- dds_vst[top_DEG_names,])

  
fname <- "FiguresX/contig_pheatmap.pdf"

pheatmap(assay(data),
           clustering_distance_cols=sampleDists,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           fontsize_row = 15,
           fontsize_col = 15,
           angle_col = 270,
           cex = 1.05,
           width = 12,
           height =12, 
          labels_col = c("Control 1","Control 2", "Control 3", "Heat stressed 1", "Heat stressed 2", "Heat stressed 3"),
           filename = fname)


sampleDistMatrix <- as.matrix(sampleDists)
sampleDists <- dist(t(assay(dds_vst)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         filename = "FiguresX/sample_pheatmap.pdf" )
  
# #################################### LRT #######################################
# 
# colData(dds_vst)
# plot_MA(res, plot_fname = "FiguresX/MAplot_LRT.pdf")
# px1 <-plot_volcano_p(res,  plot_fname = "FiguresX/EnhancedVolcano_LRT.tiff")
# px2 <- plot_volcano_padj(res, plot_fname = "FiguresX/EnhancedVolcano_padjLRT.tiff")
# 
# 
# ################################################################################
# library(DESeq2)
# library(apeglm)
# 
# library(ggpubr)
# library(devtools)
# library(ggbiplot)
# library(ggrepel)
# library(dplyr)
# library(gridExtra)
# library(tibble)
# library(magrittr)
# salmon_quants <- Sys.glob("/home/drewx/Documents/Zostera_capensis_TrX/scripts/data.enterprise/Salmon-Quant/*.sf")
# sampleTable <- read.table("/opt/DB_REF/Metadata/Z_capensis.Transcriptome2.tsv", sep = "\t", header = T)
# names(salmon_quants) <- sampleTable$Rep_name
