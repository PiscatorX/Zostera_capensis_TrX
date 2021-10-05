library(gridExtra)
library(tidyverse)
library(docstring)
library(xfun)



source("/home/drewx/Documents/Zostera_capensis_TrX/scripts/Xpression_funcs.R")

#sel_col  <-c("Sample", "FastQC_percent_duplicates", "FastQC_percent_gc","FastQC_avg_sequence_length", "FastQC_total_sequences", "FastQC_percent_fails")
#Total Raw reads
fastqc_data <- summarise_multiqc("data/multiqc/*data", "FastQC_total_sequences")
dput(colnames(fastqc_data))
(fastqc_data <- fastqc_data[c("Sample", "RawReads", "Trimmomatic", "SortmeRNA")])
(rawread_total <- sum(fastqc_data$RawReads))
(rawread_len_avg <- NA) 


#Sortmerna reads
(filtered_total <- sum(fastqc_data$SortmeRNA))
filtered_read_len <- read.table("data/infoseq/merged_infoseq_len.tsv", sep = "\t", header = T, stringsAsFactors = F)
(filtered_readlen_avg <- mean(filtered_read_len$Length))
(filtered_total <- nrow(filtered_read_len))
(filtered_readlen_min <- min(filtered_read_len$Length))
(filtered_readln_max <- max(filtered_read_len$Length))

#Assembly stats
trinity_stats <- read.table("data/Assembly.stats", sep = ":", comment.char = '#' , stringsAsFactors = F, strip.white = T) %>%
                `colnames<-`(c("Var", "value"))


#BAM mapping stats
(samtools_flagstat <- summarise_bam_metrics("data/bam_flagstats","mapped"))
(total_mapped <- sum(samtools_flagstat$mapped))

fastqc_data$Sample <- str_split_fixed(fastqc_data$Sample, "_r", n = 2)[,1]

(mapped_reads <- samtools_flagstat %>%
                rownames_to_column(var = "Sample"))

mapped_reads$Sample <- str_split_fixed(mapped_reads$Sample, "_r", n = 2)[,1]

(read_filtering <- merge(fastqc_data, mapped_reads, by = "Sample"))


{Feature <- c("Total sequencing reads (SE)",
"Read length avg. (bp)",
"Total quality filtered reads",
"Quality filtered read length avg. (bp)",
"Total Transcripts",
"Total genes",
"Transcript contig length avg.",
"Transcript N50 (bp)",
"Transcript Nx50",
"Percent GC",
"Total SE reads mapping to assembly",
"% reads mapping to assembly")}


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
total_mapped,
(total_mapped/filtered_total) * 100)

(Data_Table =  data.frame(Feature = Feature, Values = Values))

  

# filtered_data %>%
#     group_by(USA) %>%
#     summarise(mean_len =  mean(Length), max_len = max(Length), min_len = min(Length))


