#'miscellaneous functions that I use in data analysis


#' Title
#'
#' @param multiqc_path 
#' @param sel_col 
#'
#' @return
#' @export
#'
#' @examples
summarise_multiqc <- function(multiqc_path, sel_col){
  
  #filters multiqc general stats files 
  multiqc_dirs  <- Sys.glob(multiqc_path)
  
  if (length(sel_col) != 1 ){
    stop("Only one filter column required")
  }
  
  stats <- data.frame()
  
  for (sub_dir in multiqc_dirs){
    
    analysis <- basename(str_split_fixed(sub_dir ,"_", n = 2)[1])
    print(analysis)
    fname <- file.path(sub_dir,"multiqc_general_stats.txt")
    tmp_df <- read.table(fname, sep = "\t", header = T, stringsAsFactors = F) %>%
      select(c("Sample"), !!sel_col)
    tmp_df$Sample <- str_replace(tmp_df$Sample,"trim_","")
    colnames(tmp_df)[2] <- analysis
    if (nrow(stats) != 0 ){
      stats <- merge(tmp_df, stats)
    }
    else{
      stats <- data.frame(tmp_df)
    }
    
  } 
  
  return(stats)
  
}



get_var <- function(var, data =trinity_stats){
  
  trinity_stats %>%
    filter(Var == {{var}}) %>%
    pull(value)
}



#Bam file flagstat metrics
summarise_bam_metrics <- function(flagstat_dir, var){
  
  flagstat_files  <-  list.files(flagstat_dir)
  if (length(var) != 1 ){
    stop("Only one filter row should be provided")
  }
  
  flagstats <- data.frame()
  for (flagstat_fname in flagstat_files){
    tmp_df <- read.table(file.path(flagstat_dir, flagstat_fname), sep = "\t") %>%
      select(V1, V3) %>% 
      t() %>% 
      data.frame(check.names  = F, stringsAsFactors = F )
    colnames(tmp_df) <- unname(unlist(tmp_df[c("V3"),]))
    tmp_df <- tmp_df %>%  filter(row.names(tmp_df)  %in%  c("V1"))
    row.names(tmp_df) <- str_split(flagstat_fname, ".", n = 1)
    
    if (nrow(flagstats) != 0 ){
      flagstats <- rbind(flagstats, tmp_df) 
    }
    else{
      
      flagstats <- data.frame(tmp_df, check.names  = F, stringsAsFactors = F )
    }
  } 
  
  myflagstat <- flagstats %>%
    select({{var}})
  
  myflagstat[var] <- as.numeric(myflagstat[[var]])
  
  return(myflagstat)
}


# plot_dispersion <- function(dds, plot_fname ="DispEsts.tiff"){
  
#   if(interactive()){
    
#   plotDispEsts(dds)
    
#   }  
  
#   tiff(plot_fname,
#        width = 180,
#        heigh= 180,
#        units = "mm",
#        res = 1000,
#        compression = "lzw")
  

#   plotDispEsts(dds)
  
#   dev.off()
  
# }




plot_MA<- function(dds_results, plot_fname="MAplot.tiff"){
  
   if(interactive()){
    
    plotMA(dds_results)
       }  
  
   pdf(plot_fname,
        width = 180,
        heigh= 180)
  
   plotMA(dds_results)
  
   dev.off()
  
}






plot_volcano_p <- function(dds_results, plot_fname ="EnhancedVolcano.tiff", x_xlim = -5, y_xlim = 8){
  
  p <-  EnhancedVolcano(dds_results,
                        lab = rownames(dds_results),
                        x = 'log2FoldChange',
                        y = 'pvalue',
                        xlim = c(x_xlim, y_xlim))
  
  if(interactive()){
    
    print(p) 
    
  }  
  
  tiff(plot_fname,
       width = 180,
       heigh= 180,
       units = "mm",
       res = 1000,
       compression = "lzw")
  
  print(p) 
  
  dev.off()
  
 return(p)  
}

plot_volcano_padj <- function(dds_results, plot_fname ="EnhancedVolcano_padj.tiff"){
  

  p <- EnhancedVolcano(dds_results,
                       lab = rownames(dds_results),
                       x = 'log2FoldChange',
                       y = 'padj',
                       xlim=c(-6,6),
                       xlab = bquote(~Log[2]~ 'fold change'),
                       ylab = bquote(~-Log[10]~adjusted~italic(P)),
                       pCutoff = 0.001,
                       FCcutoff = 0,
                       transcriptLabSize = 4.0,
                       colAlpha = 1,
                       legend=c('NS','Log2 FC','Adjusted p-value',
                                'Adjusted p-value & Log2 FC'),
                       legendPosition = 'bottom',
                       legendLabSize = 10,
                       legendIconSize = 3.0)
  
  if(interactive()){
    
    print(p) 
    
  }  
  
  
  
  tiff(plot_fname,
       width = 180,
       heigh= 180,
       units = "mm",
       res = 1000,
       compression = "lzw")
  
  print(p) 
  
  dev.off()
  
 return(p) 

}

# subset_dds  <- function(dds_results, variable="pvalue", value=0.05, tsv_fname="resOrdered.tsv"){
  
#   dds_results <- dds_results[dds_results[[variable]] < value,]
#   resOrdered <- dds_results[order(dds_results[[variable]]),]
#   df_resOrdered  <- data.frame(resOrdered)
#   write.table(df_resOrdered, tsv_fname, sep = "\t")
  
#   return(df_resOrdered_padj_sig)
# }



# deseq_report <- function(dds, pvalueCutoff=0.05 ){
#   pvalueCutoff=0.05   
#   library(ReportingTools)
  
#   rownames(dds)  <- gsub("[^[:alnum:] ]", "_", rownames(dds))
  
#   des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2',
#                            title = 'RNA-seq analysis of differential expression using DESeq2',
#                            reportDirectory = "./reports")
  
#   publish(dds, des2Report,
#           pvalueCutoff=pvalueCutoff,
#           annotation.db="SAMSA2-nf",
#           factor = colData(dds)$condition,
#           reportDir="./reports")
  
#   finish(des2Report)
  
#   dev.off()
# }




# seasurface_ctd  <- function(ctd_data, column="day", depth_colname="Depth", depth=0.5, col_range=c(11:15), get_max=FALSE) {
  
# i = 0
# for (val in col_range){
  
#   data <- ctd_data[ctd_data[[column]]== val,]
#   data <- data[data[[depth_colname]] > depth,]
  
#   if(get_max){
#     entry_min <- min(data[[depth_colname]])
#     data <- data[data[[depth_colname]] == entry_min,]
#     }
#   else{
#      entry_min <- min(data[[depth_colname]])
#      entry_min <- sapply(entry_min, as.numeric)
#      data <- data[data[[depth_colname]] == entry_min,]
     
#      }

#   if (i == 0){
#     ssf_data  <- data.frame(matrix(data, nrow = 1))
#     colnames(ssf_data) <- colnames(data)
    
#   }
#   if  (i > 0){ 
#     ssf_data <- rbind(ssf_data,  as.numeric(data))
#   }
#   i = i + 1
# }

# ssf_data <- data.frame(sapply(ssf_data, as.numeric))
# rownames(ssf_data) <- ssf_data[[column]]

# return(ssf_data)

# }




#' @param my_dataframe Dataframe with samples(row) x measurements(columns)
#' @return mypca object list of `pca` and `var_perc` 

do_mypca  <- function(my_dataframe, scale = TRUE){



print(my_dataframe)

pca <-  prcomp(my_dataframe,  scale = scale)
print(pca)
pca_var <- pca$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var)*100,1)
print(pca$x)
pca_data <- data.frame(Sample=rownames(pca$x),
                        PC1=pca$x[,1],
                        PC2=pca$x[,2])
PC <- paste0("PC", c(1: length(pca_var_perc)))
var_perc <- data.frame(PC, pca_var_perc)
print(var_perc)
 
mypca <- list()
mypca$pca <- pca
mypca$var_perc <- var_perc
mypca$pca_data <- pca_data
return(mypca)

}



# #' @param var_df  dataframe of  
# #'
# do_mypca_varplot <- function(var_df, fname= "scree_plot.tiff" ){

# p <-ggbarplot(var_df,
#           x="PC",
#           y="pca_var_perc",
#           pallete="simpsons",
#           xlab="Principal Components",
#           title="Scree Plot",
#           ylab="Variation (%)",
#           label = round(var_df$pca_var_perc),
#           legend.title = "PC Percent variation")
# if(interactive()){
# print(p)  
  
# }

# ggsave(fname,
#        width = 180,
#        heigh= 180,
#        units = "mm",
#        dpi = 1000,
#        compression = "lzw")

#  return(p)
# }


# do_mypca_biplot <- function(pca, rowlabels="x", fname="biplot.tiff" ){

# p <-  ggbiplot(pca,labels=rownames(pca[["x"]]),
#            ellipse = TRUE,
#            labels.size = 6,
#            label.repel = FALSE,
#            circle = FALSE,
#            varname.size=5,
#            varname.adjust=1.25) +
#     theme_bw() +
#     theme(legend.direction = 'horizontal', 
#           legend.position = 'top')

# # p + scale_x_continuous(name = "Speed of cars", limits = c(0, 30)) +
# #   scale_y_continuous(name = "Stopping distance", limits = c(0, 150))

#   if(interactive()){
#     print(p)  
#     }
  
#   ggsave(fname,
#          width = 180,
#          heigh= 180,
#          units = "mm",
#          dpi = 1000,
#          compression = "lzw")
  
#   return(p)
# }



mypca_plot<- function(mypca, 
                      pca_data_index="pca_data", 
                      name = "Sample",
                      label="Sample", 
                      var_perc_index="var_perc",
                      plot_fname="PCAplot.tiff"){
 
pca_data <- mypca[[pca_data_index]]
var_perc <- mypca[[var_perc_index]]
my_colors <- colorRampPalette(c("red","blue"))(length(levels(mypca$pca_data[[label]])))

   
p <- ggplot(data=mypca$pca_data, aes_string(x="PC1", y="PC2", label=label, color = label)) +
            geom_label_repel(show.legend = FALSE, size=4)  +
            geom_point(size=1) +
            xlab(paste("PC1: ", var_perc[1], "% variance", sep="")) +
            ylab(paste("PC2: ", var_perc[2], "% variance", sep="")) +
            theme_bw() +
            theme(axis.text=element_text(size=16),
                  axis.title=element_text(size=18, face = "bold"),
                  strip.text.x = element_text(size = 14, face = "bold"),
                  axis.text.x = element_text(size = 14),
                  legend.text = element_text(size = 12),
                  legend.title = element_text(size = 14, face = "bold"),
                  panel.grid = element_blank(),
                  panel.grid.major = element_line(size=0.1)) + 
            scale_color_manual(values=my_colors, name=name)

  if(interactive()){
    print(p)  
  }
  

ggsave(plot_fname,
       width = 180,
       heigh= 180,
       units = "mm",
       dpi = 1000,
       compression = "lzw")

return(p)

}

# #' @param StHelena_Bay_CTD ctd dataframe
# #' @param depth depth at which to extract data
# #' @parm chl_max boolean whether use the chl-a max

# subsuface_chl <- function(StHelena_Bay_CTD, depth, chl_max = FALSE){
  
#   for (i in seq(1:length(days))){
#     data <- StHelena_Bay_CTD[StHelena_Bay_CTD$day == days[i],]
    
#     ordered <- data[order(data$Depth, decreasing=FALSE),]
#     row_df <- data.frame(ordered[ordered$Depth >= depth,][1,])
#     if (chl_max == TRUE ){
#       max_chl <- max(data["Chl.a"])
#       row_df$Chl.a <- max_chl
#     }
#     if( i == 1){
#       subsurf <- row_df
#     }
#     else{
#       subsurf <-rbind(subsurf, row_df)
#     }
#   }
#   return(subsurf)
# }



# merge_data <- function(source_dir, var_col=3, header = FALSE, val_col=1, sep  = "\t", ext = "tsv"){
  
#   files =  Sys.glob(file.path(source_dir, paste0("*",ext)))
#   i = 1
#   for (datafile in files){
#     data_df = read.delim(datafile, stringsAsFactors = FALSE, sep = sep, header = header)[, c(var_col, val_col)]
#     name <- basename(datafile)
#     var <- paste0("data.",i)
#     colnames(data_df) <- c(var, name )
    
#     if (i == 1){
#       merged_df  <- data_df
#       merge_col = var
#     }
#     else{
#       merged_df <- merge(merged_df, data_df, by.x = merge_col,  by.y = var, all = TRUE)
#     }
#     i = i + 1
#   }
#   merged_df[is.na(merged_df)] <- 0
#   return(merged_df)
# } 

