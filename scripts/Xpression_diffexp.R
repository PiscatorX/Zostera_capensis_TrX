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
library(xtable)


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


############ significant DE analysis ###########################################
#https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
res_dat <- res %>%
           data.frame %>%
           rownames_to_column(var="SwissProt")

write.table(res_dat, file = "data.enterprise/Blast_top_hits/TrX_padj1.tsv", row.names = FALSE, quote = FALSE, sep ="\t", col.names = T)

  
signf_genes <- function(padj_cutoff){
  
    res_sig <- res%>%
               data.frame() %>%
               rownames_to_column(var="ACC") %>%
               filter(padj < {padj_cutoff}) %>%
               filter(abs(log2FoldChange) > 1) %>%
               arrange(desc(log2FoldChange)) %>%
               as_tibble

return(res_sig)

}

(sig_genes1 <- signf_genes(padj_cutoff =  0.05))
colnames(sig_genes1)
colnames(uniprot_ref)
dim(sig_genes1)
(top_DEG_names <- sig_genes1 %>% pull(ACC))

DEG_data <- merge(sig_genes1, uniprot_ref, by = "ACC", all.x = T) %>%
            select(ACC, name, log2FoldChange, lfcSE, padj) %>%
            arrange(desc(log2FoldChange))
            
dim(DEG_data)

print(xtable(DEG_data, digits= 12, type = "latex"), file = "SuppInfoX/DEG_data.tex", include.rownames = F)

#top_DEG_order <- order(rowMeans(counts(dds_dexp[top_DEG_names,], normalized=TRUE)), decreasing=TRUE)
  
############################ Heatmap  ##########################################
  
colData(dds_vst)
metadata_col <- as.data.frame(colData(dds_vst)[,c("Rep_name")])
colnames(metadata_col) <-  c("Replicate")
metadata_row <- as.data.frame(colData(dds_vst)[,c("Rep_name")])
colnames(metadata_row) <-  c("Sample")
rownames(metadata_row) <- sampleTable$Readname
  
colnames(dds_vst) <- sampleTable$Rep_name
(data <- dds_vst[top_DEG_names,])

  
fname <- "FiguresX/contig_pheatmap.pdf"

heatmapdata <- assay(data) %>%
               data.frame() %>%
               rownames_to_column(var = "ACC")



dim(heatmapdata)
colnames(uniprot_ref)
heatmapdata_ref <- merge(heatmapdata, uniprot_ref, by = "ACC") 

heatmapdata_ref %>%
  group_by(name) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
  
make.unique(as.character(heatmapdata_ref$name))


heatmap  <- heatmapdata_ref %>%
            select(-c(ACC, name))

dim(heatmapdata_ref)

summary(rowMaxs(as.matrix(heatmap)))

rownames(heatmap) <- make.unique(as.character(heatmapdata_ref$name))

sampleDistMatrix <- as.matrix(sampleDists)
sampleDists <- dist(t(assay(dds_vst)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         filename = "FiguresX/sample_pheatmap.pdf" )


pheatmap(heatmap,
           clustering_distance_cols=sampleDists,
           cluster_cols = F,
           cluster_rows = T,
           show_rownames = F,
           fontsize_row = 15,
           cellwidth = 25,
           cellheight = 5,
           fontsize_col = 15,
           angle_col = 90,
           treeheight_row = 75,
           cex = 1.1,
           legend_breaks = seq(6,14,1), 
           legend_labels = seq(6,14,1), 
           labels_col = c("Control 1","Control 2", "Control 3", "Heat stressed 1", "Heat stressed 2", "Heat stressed 3"),
           filename = fname)


  
