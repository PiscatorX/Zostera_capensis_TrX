library(reticulate)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(ti)
library(qpcR)
library(stringi)
library(gridExtra)
library(cowplot)
library(reshape2)


fastacounter  <- import('fastacounta')

f_list <- Sys.glob("/home/drewx/Documents/Zostera_capensis_TrX/scripts/data.enterprise/Assemblies/*.fasta")

len_stats <- do.call(qpcR:::cbind.na, lapply(f_list, function(fname){
  
  return(fastacounter$len_stats(fname, "fasta")[2])
}))


filenames <-  data.frame(fnames = cbind(lapply(f_list,  basename))) %>%
             separate(fnames, sep = '\\.', c("algorithm", "extension"))

names(len_stats) <- filenames$algorithm

################### sequence statistics ##########################################
len_data <-  do.call(rbind, lapply(len_stats, function(col){return(sum(!is.na(col)))}))  %>%
               data.frame() %>%
              `colnames<-`("transcripts") %>%
               rownames_to_column(var  = "Assembler") %>%
               arrange(transcripts)

len_data$Assembler <- factor(len_data$Assembler, levels = len_data$Assembler)

################################################################################

options(scipen=999)

ggplot(len_data, aes(x = Assembler, y = transcripts)) +
          geom_bar(stat = "identity") +
          ylim(limits = c(0,600000 )) +
          ylab("Total transcripts") +
          xlab("Assembler algorithm") +
            theme(panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_text(size=18, colour = "black"),
                  axis.title = element_text(size=18, colour = "black"),
                  panel.border = element_rect(colour = "black", fill=NA, size=1))



ggsave("SuppInfoX/assemblers.pdf")      

################################################################################

melt_len_stats <-  melt(len_stats)

p1 <- ggplot(melt_len_stats %>% filter(variable != "EvidentialGene") , aes(x = value)) + 
            geom_bar() + 
            ylab("Count") +
            xlab("Length (bp)") +
            theme(panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_text(size=14, colour = "black"),
                  axis.title = element_text(size=14, colour = "black"),
                  strip.background = element_rect(colour = "black", fill=NA, size=1),
                  panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                  geom_rug() +
            facet_wrap(~variable, scales = "free_y", ncol = 2)
      


p2 <- ggplot(melt_len_stats %>% filter(variable == "EvidentialGene") , aes(x = value)) + 
        geom_bar() + 
        ylab("Count") +
        xlab("Length (bp)") +
        theme(panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text = element_text(size=14, colour = "black"),
              axis.title = element_text(size=14, colour = "black"),
              strip.background = element_rect(colour = "black", fill=NA, size=1),
              panel.border = element_rect(colour = "black", fill=NA, size=1)) +
              geom_rug()
  



gt <- arrangeGrob(p1,                               
                  p2, 
                  ncol = 1, nrow = 2)

p <- as_ggplot(gt) +                      
     draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0), y = c(1, 0.5))

p

ggsave("SuppInfoX/Evigen_len_dist.pdf")      

