library(gridExtra)
library(tidyverse)
library(docstring)
library(xfun)



source("/home/drewx/Documents/Zostera_capensis_TrX/scripts/Xpression_funcs.R")

#sel_col  <-c("Sample", "FastQC_percent_duplicates", "FastQC_percent_gc","FastQC_avg_sequence_length", "FastQC_total_sequences", "FastQC_percent_fails")

fastqc_data <- summarise_multiqc("data/multiqc/*data", "FastQC_total_sequences")
dput(colnames(fastqc_data))
fastqc_data <- fastqc_data[c("Sample", "RawReads", "Trimmomatic", "SortmeRNA")]
(Total_reads <- sum(fastqc_data$RawReads))

#Raw reads
rawread_data <- read.table("data/infoseq/merged_infoseq.tsv", sep = "\t", header = T, stringsAsFactors = F)
(rawread_len_avg <- mean(rawread_data$Length) %>% round(0) ) 
(rawread_total <- nrow(rawread_data))


#Sortmerna reads
filtered_data <- read.table("data/infoseq/merged_infoseq.tsv", sep = "\t", header = T, stringsAsFactors = F)
filtered_data$USA  <- as.data.frame(str_split_fixed(filtered_data$USA,":", n = 4))[3]
(filtered_total <- nrow(filtered_data))
(filtered_readlen_avg <- mean(filtered_data$Length) %>% round(0))
(filtered_readlen_min <- min(filtered_data$Length))
(filtered_readln_max <- max(filtered_data$Length))

#Assembly stats
trinity_stats <- read.table("data/Assembly.stats", sep = ":", comment.char = '#' , stringsAsFactors = F, strip.white = T) %>%
                `colnames<-`(c("Var", "value"))


#Samtools Mapping stats

list.files(file.path(bam_dir, "*.metrics"))


summarise_bam_metrics <- function(bam_dir, var_col){
  
  metric_files  <-  list.files(bam_dir)
  print(metric_files)
  
  if (length(var_col) != 1 ){
    
    stop("Only one filter column required")
  
  }

  metrics <- data.frame()

  for (read_metric in metric_files){
    
      tmp_df <- read.table(fname, sep = "\t", header = T, stringsAsFactors = F) %>%
       select(c("Sample"), !!sel_col)
  #   tmp_df$Sample <- str_replace(tmp_df$Sample,"trim_","")
  #   colnames(tmp_df)[2] <- analysis
  #   if (nrow(stats) != 0 ){
  #     stats <- merge(tmp_df, stats) 
  #   }
  #   else{
  #     stats <- data.frame(tmp_df)
  #   }
  #   
  # } 
  # 
  # return(stats)
  
}


summarise_bam_metrics("/home/drewx/Documents/Zostera_capensis_TrX/scripts/data/samtools-metrics", "")



read.table(file.path( IonXpressRNA_004_rawlib.metrics" )


Feature <- c("Total sequencing reads (SE)",
"Read length avg. (bp)",
"Total quality filtered reads",
"Quality filtered read length avg. (bp)",
"Total Transcripts",
"Total genes",
"Transcript contig length avg.",
"Transcript N50 (bp)",
"Transcript Nx50",
"Percent GC",
"Total SE reads mapping to assembly")


Values = c(rawread_total,
rawread_len_avg,
filtered_total,
filtered_readlen_avg,
get_var("Total trinity transcripts"),
get_var("Total trinity genes"),
get_var("Average contig")[1],
get_var("Contig N50")[1],
NA,
get_var("Percent GC"),
NA)

(Data_Table =  data.frame(Feature = Feature, Values = Values))

  




filtered_data %>%
    group_by(USA) %>%
    summarise(mean_len =  mean(Length), max_len = max(Length), min_len = min(Length))


