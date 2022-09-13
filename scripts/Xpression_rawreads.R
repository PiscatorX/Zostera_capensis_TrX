library(gridExtra)
library(tidyverse)
library(docstring)
library(xfun)
library(xtable)


source("/home/drewx/Documents/Zostera_capensis_TrX/scripts/Xpression_funcs.R")

#sel_col  <-c("Sample", "FastQC_percent_duplicates", "FastQC_percent_gc","FastQC_avg_sequence_length", "FastQC_total_sequences", "FastQC_percent_fails")
#Total Raw reads
fastqc_data <- summarise_multiqc("data/multiqc/*data", "FastQC_total_sequences")
dput(colnames(fastqc_data))
(fastqc_data <- fastqc_data[c("Sample", "RawReads", "Trimmomatic", "SortmeRNA")])
(rawread_total <- sum(fastqc_data$RawReads))


#Raw length distribution
#rawread_len <- read.table("data.enterprise/infoseq/len.dist", sep = "\t", header = T, stringsAsFactors = F)
#save(rawread_len,file="data.enterprise/infoseq/rawread_len.Rda")
load("data.enterprise/infoseq/rawread_len.Rda")
(rawread_len_avg <- mean(rawread_len$Length))
summary(rawread_len$Length)

rawdata_len <- summarise_multiqc("data/multiqc/*data", "FastQC_avg_sequence_length")
summary(rawdata_len$RawReads)

sum(rawdata_len$RawReads)

summary(rawdata_len$SortmeRNA)

(readsums <- colSums(fastqc_data[c("RawReads","Trimmomatic","SortmeRNA")]))

(readsums/67780422)

#Sortmerna reads
(filtered_total <- sum(fastqc_data$SortmeRNA))
#filtered_read_len <- read.table("data/infoseq/merged_infoseq.tsv", sep = "\t", header = T, stringsAsFactors = F)
#save(filtered_read_len,file="data/infoseq/filtered_read_len.Rda")
load("data/infoseq/filtered_read_len.Rda")
read_len <- filtered_read_len %>%
            group_by(Length) %>% 
            summarise(Count = n())

summary(filtered_read_len)
(filtered_readlen_avg <- mean(filtered_read_len$Length))
(filtered_total <- nrow(filtered_read_len))
(filtered_readlen_min <- min(filtered_read_len$Length))
(filtered_readln_max <- max(filtered_read_len$Length))

#Assembly stats
evigene_stats <- read.table("data.enterprise/Evigene_alt.assembly_stats", sep = ":", comment.char = '#' , stringsAsFactors = F, strip.white = T) %>%
                `colnames<-`(c("Var", "value"))


#BAM mapping stats
(samtools_flagstatx <- summarise_bam_metrics("data.enterprise//bam_flagstats","mapped"))
(total_mappedx <- sum(samtools_flagstatx$mapped))

fastqc_data$Sample <- str_split_fixed(fastqc_data$Sample, "_r", n = 2)[,1]

(mapped_readsx <- samtools_flagstatx %>%
                rownames_to_column(var = "Sample"))

mapped_readsx$Sample <- str_split_fixed(mapped_readsx$Sample, "_r", n = 2)[,1]

(read_filteringx <- merge(fastqc_data, mapped_readsx, by = "Sample"))

100 *  round(sum(read_filteringx$mapped)/sum(read_filteringx$SortmeRNA),1)


(read_filteringx$mapping_rate <-  round(100*(read_filteringx$mapped/read_filteringx$SortmeRNA),1))

print(xtable(read_filteringx, digits=7, type = "latex"), file = "SuppInfoX/read_filtering.tex", include.rownames = F)


################################################################################

data <-  read.table("SuppInfoX/fastqc_sequence_length_distribution_plot.tsv", sep = "\t", header = T, check.names = F)

cnames <- colnames(data)

colnames(data) <-NULL

data.frame(t(data))

