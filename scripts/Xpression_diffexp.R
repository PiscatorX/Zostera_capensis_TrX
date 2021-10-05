#!/usr/bin/env Rscript
rm(list=ls())
cat("\f")

# library(DESeq2)
# library(tximport)
# library(BiocParallel)
# library(DESeq2)
# library(apeglm)
# library(EnhancedVolcano)
# library(ggpubr)
# library(devtools)
# library(ggbiplot)
# library(ggrepel)
# library(dplyr)
# library(gridExtra)
# library(tibble)
# library(magrittr)
# library(pheatmap)


source("/home/drewx/Documents/Zostera_capensis_TrX/scripts/Xpression_funcs.R")


trX <- Sys.glob("/home/drewx/Documents/Zostera_capensis_TrX/scripts/data/Salmon-Quant/*.sf")

sampleTable <- data.frame()
# condition_d2  <- data.frame(c("D1", "D1", "D1", "D2"))
# colnames(condition_d2)  <- c("DaysSinceSamplingStart")
sampleTable <- data.frame(condition = condition_d2


# #tx_d2 <- tximport(files_d2, type="salmon", txOut=T)

# tx_tibble <- tx_d2$counts  %>% data.frame()  %>%  cbind(row.names(tx_d2$counts))  %>% as_tibble()
# samplenames  <- c(paste0("D1.rep", seq(3)), "D3", "contig_id")
# colnames(tx_tibble) <-  samplenames

# counts <- tx_tibble[order(-tx_tibble$D1.rep1),]
# write.table(counts, file = "transcipt_counts.tsv", row.names = FALSE, quote = FALSE, sep ="\t")

# rownames(sampleTable) <- colnames(tx_d2$counts)
# dds_dataset <- DESeqDataSetFromTximport(tx_d2, sampleTable, ~DaysSinceSamplingStart)
# dds_LRT <- DESeq(dds_dataset, minReplicatesForReplace=Inf)
# res_LRT = results(dds_LRT)

# #################################################### LRT #############################################################
# dds_vst_LRT <- vst(dds_LRT, blind=FALSE)
# colData(dds_vst_LRT)

# metadata_col <- as.data.frame(colData(dds_LRT)[,c("DaysSinceSamplingStart")])
# colnames(metadata_col) <-  c("Sampling Day")
# metadata_row <- as.data.frame(colData(dds_LRT)["DaysSinceSamplingStart"])
# colnames(metadata_row) <-  c("Sample")
# rownames(metadata_row)  <-  samplenames[-c(5)]

# # plot_dispersion(dds_LRT, plot_fname = "DispEsts_LRT.tiff")
# # plot_MA(res_LRT, plot_fname = "MAplot_LRT.tiff")
# # px1 <-plot_volcano_p(res_LRT,  plot_fname = "EnhancedVolcano_LRT.tiff")
# px2 <- plot_volcano_padj(res_LRT, plot_fname = "EnhancedVolcano_padjLRT.tiff")


# ########################################################### significant DE analysis ###########################################
# #https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
  
# res_LRT_sig <- res_LRT%>%
#   data.frame() %>%
#   rownames_to_column(var="RefSeq") %>% 
#   filter(padj < 0.001) %>%
#   filter(abs(log2FoldChange) > 2) %>%
#   arrange(desc(log2FoldChange)) %>%
#   as_tibble

# write.table(res_LRT_sig, file = "SHB_MTx_padj.tsv", row.names = FALSE, quote = FALSE, sep ="\t")






# top_DEG_data <- res_LRT_sig %>% top_n(50, baseMean) 
# top_DEG_names <- top_DEG_data %>% pull(RefSeq)
# top_DEG_order <- order(rowMeans(counts(dds_LRT[top_DEG_names,], normalized=TRUE)), decreasing=TRUE)
  
#   ########################################################### heatmap  ##################################################
  
#   dds_vst_LRT <- vst(dds_LRT, blind=FALSE)
#   colData(dds_vst_LRT)
  
#   metadata_col <- as.data.frame(colData(dds_LRT)[,c("DaysSinceSamplingStart")])
#   colnames(metadata_col) <-  c("Sampling Day")
#   metadata_row <- as.data.frame(colData(dds_LRT)["DaysSinceSamplingStart"])
#   colnames(metadata_row) <-  c("Sample")
#   rownames(metadata_row) <- samplenames[-c(5)]
  
  
#   data <- dds_vst_LRT[top_DEG_names,]
  
#   fname <- "contig_pheatmap.pdf"
  
#   normalised_counts <- assay(data)[top_DEG_order,]
  
#   write.table(normalised_counts, 
#               "shb_normalised.tsv",
#               quote = FALSE,
#               sep = "\t")
  
#   pheatmap(assay(data)[top_DEG_order,],
#            cluster_cols = FALSE,
#            cluster_rows = FALSE,
#            angle_col = 0,
#            cex = 1.05,
#            filename = fname)
  
  
#   ##################################################### PCA analysis ############################################################## 
#   mypca_LRT <- list()
#   mypca_LRT$pca_data <- plotPCA(dds_vst_LRT, intgroup ="DaysSinceSamplingStart", returnData=TRUE)
#   mypca_LRT$var_perc <- round(100 * attr(mypca_LRT$pca_data, "percentVar"))
#   mypca_LRT$pca_data$name <- factor(c("D1", "D1", "D1", "D3"))
#   px3 <- mypca_plot(mypca_LRT, label="name", plot_fname = "PCAplot_LRT.tiff")
  
