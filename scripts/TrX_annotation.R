library(RColorBrewer)
library(treemapify)
library(reticulate)
library(tidyverse)
library(reshape2)
library(magrittr)
library(ggbreak) 
library(cowplot)
library(ggrepel)
library(xtable)
library(scales)
library(ggpubr)
library(ggsci)
library(aplot)
library(grid)


####################################### Swiss-Prot ###################################
fastacounter  <- import('fastacounta')
transcripts <- fastacounter$len_stats("data.enterprise/Assemblies/EvidentialGene.fasta", "fasta")
(total_transcripts <- length(transcripts$seq_id))
summary(transcripts$len)

Blast_Swissprot  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.ids", sep = "\t") %>%
                    `colnames<-`(c("transcript","Uniprot_ID"))

(total_Swissprot_transcripts <- length(Blast_Swissprot$transcript))

uniprot_sprot = read.csv("data.enterprise/Blast_top_hits/uniprot_sprot.ref", sep = "\t", header = F) %>%
                `colnames<-`(c("Uniprot_ID","name"))


(uniq_Swissprot <- Blast_Swissprot %>%
    group_by(Uniprot_ID) %>%
    summarise(count = n())) %>%
  arrange(desc(count))


(total_uniq_Swissprot <- length(uniq_Swissprot$Uniprot_ID))

uniq_Swissprot_ref <- merge(uniq_Swissprot, uniprot_sprot, by = "Uniprot_ID", all.x = T) %>%
  arrange(desc(count))


swissprot_plot  <- uniq_Swissprot_ref %>%
                    mutate(label = paste0(name)) %>%
                    mutate(perc = round(count/sum(count)*1000,1)) %>%
                    arrange(desc(count)) %>%
                    head(25) %>%
                    ggbarplot(x = "label", y = "perc",
                          fill = "#54b2a9",
                          color = "white",
                          palette = "jco",
                          sort.val = "asc",
                          ylab = "Annotated transcripts (%)",
                          xlab = "Swiss-Prot proteins",
                          orientation = "horiz"
                    ) + theme_bw() +
                    theme(legend.background = element_rect(colour = "darkgray"),
                      axis.text=element_text(size=12),
                      panel.grid = element_line(colour = "black",linetype="dashed",size=0.1),
                      axis.title = element_text(size=12),
                      legend.text = element_text(size=12)) + 
                    scale_y_continuous(breaks = seq(0,12.5, 2.5), labels = seq(0,12.5,2.5), limits = c(0, 12.5))


swissprot_plot
 
ggsave("SuppInfoX/Swiss_Prot.pdf")
print(xtable(swissprot_plot$data %>% select(-label), digits=7, type = "latex"), file = "SuppInfoX/SwissProt.tex", include.rownames = F) 

############################# Swiss-Prot GO #####################################
Swissprot_GO_all  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.map", sep = "\t") %>%
                     `colnames<-`(c("transcript","Uniprot_ID","GO_accession","go_annotation","source")) 

length(unique(Swissprot_GO_all$Uniprot_ID))

Swissprot_GO_raw <-  Swissprot_GO_all %>% 
                     dplyr::select(transcript, GO_accession)


dim(Swissprot_GO_raw)
Swissprot_GO_parsed  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx_map_clean.tsv", sep = "\t") %>%
                         `colnames<-`(c("transcript","Uniprot_ID","GO_prefix","GO_ID","domain","term"))
dim(Swissprot_GO_parsed)

head(Swissprot_GO_raw$transcript) ==  head(Swissprot_GO_parsed$transcript)
tail(Swissprot_GO_raw$transcript) == tail(Swissprot_GO_parsed$transcript)
length(Swissprot_GO_raw$transcript) == length(Swissprot_GO_parsed$transcript)

Swissprot_GO_parsed <- Swissprot_GO_parsed %>%
                       dplyr::select(Uniprot_ID, GO_prefix, GO_ID, domain, term)

Swissprot_GO_tmp  <- cbind(Swissprot_GO_raw, Swissprot_GO_parsed)
length(Swissprot_GO_tmp$transcript) == length(Swissprot_GO_raw$transcript)
noquote(paste(colnames(Swissprot_GO_tmp), collapse = ','))

Swissprot_GO  <- Swissprot_GO_tmp %>%  
                 dplyr::select(transcript, GO_accession, transcript, Uniprot_ID, domain, term)

GO_annotation  <- Swissprot_GO %>%
                  group_by(GO_accession) %>%
                  summarise(count = n()) %>%
                  arrange(desc(count))

(total_GO_terms <- length(GO_annotation$GO_accession))

Swissprot_transcripts <- Swissprot_GO %>% 
                          group_by(transcript) %>%
                          summarise(count = n())

(total_Swissprot_transcripts <- length(Swissprot_transcripts$transcript))

GO_transcripts <-  Swissprot_GO %>%
                    group_by(transcript) %>%
                    summarise(count = n()) %>%
                    arrange(desc(count))

(total_GO_transcripts <- length(GO_transcripts$transcript))

plant_GO_slim <-  read.csv("data.enterprise/Blast_top_hits/goslim_plant.tsv", sep = "\t",  header = F) %>%
                 `colnames<-`(c("GO_accession", "goslim_domain", "goslim_term"))

Swissprot_GO_plant <- Swissprot_GO[Swissprot_GO$GO_accession %in% plant_GO_slim$GO_accession,]  

all_transcripts <- Swissprot_GO_all %>%
                   distinct(transcript)

dim(all_transcripts)

Plant_transcripts <- Swissprot_GO_plant %>%
                     distinct(transcript)

(plant <- length(Plant_transcripts$transcript))  

`%nin%` <- Negate(`%in%`)

non_plant_transcripts <-  all_transcripts$transcript[all_transcripts$transcript %nin%  Plant_transcripts$transcript]

(no_plant <- length(non_plant_transcripts))

non_plantGO_transcripts <- length(Swissprot_GO$transcript) - length(Swissprot_GO_plant$transcript)

dim(non_plantGO_transcripts)

GOslim_annotation  <- Swissprot_GO_slim %>%
                      group_by(GO_accession) %>%
                      summarise(count = n()) %>%
                      arrange(desc(count))

sum(GOslim_annotation$count)

(total_GOslim_terms <- length(GOslim_annotation$GO_accession))


Swissprot_GO_slim  <- merge(Swissprot_GO_plant, plant_GO_slim, by = "GO_accession") 

write.table(Swissprot_GO_slim, 
            file = "data.enterprise/Blast_top_hits/Swissprot_GO_slim.tsv",
            row.names = FALSE, 
            quote = FALSE, 
            sep ="\t", 
            col.names = T)



length(unique(Swissprot_GO_slim$Uniprot_ID))


noquote(colnames(Swissprot_GO_slim))

plant_GO_slim %>%
              arrange(goslim_domain, goslim_term)

groupGO <- plant_GO_slim %>% 
           group_by(goslim_domain) %>%
           summarise(count = n())


(Swissprot_GO_slim_count <- Swissprot_GO_slim %>%
                          group_by(goslim_domain, goslim_term) %>%
                          summarise(count = n()) %>%
                          arrange(desc(count)) %>%
                          top_n(15, count))


Swissprot_GO_slim %>%
          group_by(goslim_domain) %>%
          summarise(count = n())


dim(Swissprot_GO_slim)

unique(Swissprot_GO_slim_count$goslim_domain)

(GO_sorted <-  Swissprot_GO_slim %>%
               group_by(GO_accession, goslim_domain, goslim_term) %>%
               summarise(count = n()) %>%
               arrange(desc(count)) %>%
               arrange(goslim_domain, desc(count))
)


print(xtable(GO_sorted, digits=7, type = "latex"), file = "SuppInfoX/GO_sorted.tex", include.rownames = F)            

Swissprot_GO_slim_count$goslim_domain <- factor(Swissprot_GO_slim_count$goslim_domain, levels = c("biological_process", "molecular_function","cellular_component"))

swissprotGO <- ggbarplot(Swissprot_GO_slim_count, x = "goslim_term", y = "count",
                        fill = "goslim_domain",               
                        color = "white",
                        palette = "jco",
                        sort.val = "asc",  
                        ylab = "Swiss-Prot proteins",
                        xlab = "GO term",
                        legend.title="GO category",
                        orientation = "horiz") + 
                theme_bw() +
                theme(legend.background = element_rect(colour = "darkgray"),
                      legend.title = element_text(size=16),
                      legend.position = c(0.775, 0.1),
                      axis.text.y = element_text(size=16),
                      axis.text.x =element_text(size=16, angle = 270, hjust = 0, vjust = 1),
                      panel.grid = element_line(colour = "darkgray", linetype="dashed", size=0.1),
                      axis.title = element_text(size=16),
                      legend.text = element_text(size=16)) + 
                scale_y_continuous(breaks =seq(0,30000,5000)) +
                scale_fill_discrete(labels = c("Biological process", "Molecular function","Cellular component"))

swissprotGO 
ggsave("FiguresX/go_terms.pdf", width = 30, height  = 30, units = c("cm"))

############################ Swiss-Prot Taxonomy ###############################

Swissprot_species  <- read.csv("data.enterprise/Blast_top_hits/uniprot_sprot.species_map", sep = "\t") %>%
                        `colnames<-`(c("Uniprot_ID","species"))

Swissprot_GO_parsed  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx_map_clean.tsv", sep = "\t") %>%
  `colnames<-`(c("transcript","Uniprot_ID","GO_prefix","GO_ID","domain","term"))
dim(Swissprot_GO_parsed)
colnames(Swissprot_GO_parsed)

swissprot_count  <-  Swissprot_GO_parsed %>%
                     group_by(Uniprot_ID) %>%
                     summarise(count = n()) %>%
                     arrange(desc(count)) 

dim(swissprot_count) 

                 
species_map <- merge(swissprot_count, Swissprot_species, by = "Uniprot_ID", all.x =  T) %>%
               arrange(species)
dim(species_map)

species_map[species_map$species=="Acanthamoeba castellanii",] %>% 
                                            group_by(species) %>%
                                            summarise(count = sum(count))

species_counts <- species_map %>%
                    group_by(species) %>%
                    summarise(count = sum(count)) %>%
                    mutate(perc = round(100*count/sum(count),2)) %>%
                    mutate(label = paste0(species, " [",perc,"%]")) %>%
                    arrange(desc(perc)) %>%
                    head(15)


other_perc <- 100 - sum(species_counts$perc)
other <- c("Other taxa", sum(species_map$count) - sum(species_counts$count), other_perc, paste0("Other", " [",other_perc,"%]"))
species_counts <- rbind(species_counts, other)
species_counts$perc <- as.numeric(species_counts$perc)

species_counts$species <- factor(species_counts$species, levels =  species_counts$species)

getPalette = colorRampPalette(brewer.pal(9, "Set2"))

species_counts %>%
      arrange(desc(count))

ggplot(data = species_counts, aes(area = perc, fill = label, label = species)) +
       geom_treemap() +
       scale_fill_manual(name = "Swiss-Prot species annotation", values = getPalette(length(species_counts$species))) +
  geom_treemap_text(fontface = "italic") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title =  element_text(size=14, colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 

       
ggsave("SuppInfoX/SwissProt_taxonomy.pdf")

######################################### eggNOG ################################

fixCOG  <- import('fixCOG')

eggNOG <- read.table("data.enterprise/eggNOG/Evigene.emapper.annotations", fill=TRUE, sep = "\t", header = T)
COG_category <- read.table("data.enterprise/eggNOG/fun-20.tab", fill=TRUE, sep = "\t") %>%
  `colnames<-`(c("COG_category","Hexadecimal","Category_description"))

colnames(eggNOG)
dim(eggNOG)
COG_count <- eggNOG[eggNOG$COG_category != '-',]
dim(COG_count)

COG_category_count <- eggNOG %>%
  select(query, COG_category, Preferred_name, Description) %>%
  group_by(COG_category) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

fix_COG_category_count <- fixCOG$getCOGsingletons(COG_category_count)
fix_COG_category_count$count <- unlist(fix_COG_category_count$count)

splice_COG_category_count <-  fix_COG_category_count %>% 
  group_by(COG_category) %>%
  summarise(countx = sum(count))



COG <- merge(COG_category, splice_COG_category_count)


COG$label <- paste(COG$Category_description, COG$COG_category,  sep = ": ")
COG$Hexadecimal <- paste("#", COG$Hexadecimal, sep ="")

COG$label <- factor(COG$label, levels = rev(COG$label))

COG %>%
  arrange(desc(countx)) %>%
  select(COG_category, Category_description, countx)

COG_category_count %>%
  arrange(desc(count))



COGplot <- ggplot(COG, aes(x = label, y = countx)) +
  geom_bar(stat= "identity", fill = "#00AFBB") +
  scale_fill_manual() +
  xlab("COG category") +  
  ylab("Number of transcripts Log(count)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size=16),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top =  element_blank(),
        axis.text.x = element_text(size=16, angle = 270, vjust = 0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16, hjust = 1, vjust = 0),
        legend.position = "top",
        legend.text = element_text(size=16, colour = "black"),
        legend.background = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(colour = "black", fill=NA, size=1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
        scale_y_break(c(2000, 4500)) +
  scale_y_continuous(breaks = c(seq(0,2000, 500), seq(4500, 5000,250))) +
  coord_flip() 


COGplot 

ggsave("FiguresX/eggNOG_COG.pdf")

############################### Pfam ########################################### 
pfam_results <- read.csv("data.enterprise/transdcoder/pfam.domtblout.tsv",  sep = "\t") %>%
  select(-accession.1) 

colnames(pfam_results)

(pfam_results_count  <- pfam_results %>%
    group_by(query.name) %>%
    summarise(count = n()) %>%
    arrange(desc(count)))

dim(pfam_results_count)

pfam_results_count %>%
  distinct(count) %>%
  arrange(desc(count))

pfam_results_filter <-  pfam_results %>%
  filter(E.value <= 1e-5) %>%
  group_by(query.name) %>% 
  slice(which.min(E.value)) %>%
  ungroup() %>% 
  mutate(accession = str_split_fixed(accession, '[.]',  2)[,1])


(pfam_results_filter_count <- pfam_results_filter %>%
    group_by(query.name) %>%
    summarise(count = n()))

dim(pfam_results_filter_count)

pfam_results_filter_count %>% 
  distinct(count) 

lapply(list(pfam_results, pfam_results_filter), function (df){
  dim(df)
})

head(pfam_results_filter)

pfam_results_counts <- pfam_results_filter  %>%
  group_by(accession) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

(pfam_deflines <- read.csv("data.enterprise/Pfam/Pfam.def", sep = "\t",  header = F) %>%
    `colnames<-`(c("accession" ,  "name",  "description")) %>%
    mutate(accession = str_split_fixed(accession, '[.]',  2)[,1])) %>% head()


pfam_deflines$accession[pfam_deflines$accession %in% pfam_results_filter$accession]

str(pfam_deflines$accession)
str(pfam_results$accession)

lapply(list(pfam_deflines, pfam_results_filter), function(df) {
  dim(df)
  colnames(df)  
})


pfam_results_defs <- merge(pfam_deflines, pfam_results_filter)

pfam_results_defs %>%
  group_by(query.name) %>%
  summarise(count = n())%>%
  distinct(count) 

`%nin%` <- Negate(`%in%`)

dropped <- pfam_results_count[pfam_results_count$query.name %nin% pfam_results_defs$query.name, ]

pfam_results[pfam_results$query.name %in% dropped$query.name, ] %>%
  select(target.name,  accession, E.value, query.name) %>% 
  arrange(E.value)

lapply(list(pfam_deflines, pfam_results_filter, pfam_results_defs), function(df) {
  
  print(dim(df))
  print(colnames(df))
  
  return(NULL)
})

write.table(pfam_results_defs, 
            file = "data.enterprise/Pfam/pfam_results_defs.tsv",
            row.names = FALSE, 
            quote = FALSE, 
            sep ="\t", 
            col.names = T)


colnames(pfam_results_defs)


pfam_results_defs %>% head() 


pfam_plot <-  pfam_results_defs %>%
  select(accession, name, description,  query.name) %>%
  group_by(name,description) %>%
  summarise(count = n()) %>%
  mutate(label = paste0(description,' [',name,']')) %>%
  ungroup() %>%
  mutate(perc = round(count/sum(count)*1000,2)) %>%
  arrange(desc(count)) %>%
  head(25) %>%
  ggbarplot(x = "label", y = "perc",
            fill = "#0392cf",
            color = "white",
            palette = "jco",
            sort.val = "asc",
            ylab = "Transcripts (%)",
            xlab = "Pfam domains",
            legend.title="GO Domain",
            orientation = "horiz"
  ) + theme_bw() +
  theme(legend.background = element_rect(colour = "darkgray"),
        axis.text=element_text(size=16),
        panel.grid = element_line(colour = "black",linetype="dashed",size=0.1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16)) + 
  scale_y_continuous(breaks =seq(0,30,5))



pfam_plot

############################## panel plot Transcripts ##########################

design <- "
   12
   13
  "

plot_list(gglist=list(swissprotGO, COGplot, pfam_plot), design=design)

ggsave("FiguresX/TransciptAnnotation.pdf", width = 24, height = 14)

######################################### KEGG #################################
ko_ids <- read.csv2("data.enterprise/KAAS/ko.ids", sep = "\t", strip.white = T, header = F) %>%
       `colnames<-`(c("KO"))

dim(ko_ids)

ko_count <- ko_ids %>% 
            group_by(KO) %>%
            summarise(count = n()) %>%
            arrange(desc(count))
    
dim(ko_count)
ko_map <- read.csv("data.enterprise/KAAS/kegg.name_category", sep = "\t", strip.white = T, header = F) %>%
          `colnames<-`(c("KO","name", "category")) %>%
           distinct()
dim(ko_map)

ko_map %>% 
          group_by(KO) %>%
          summarise(count = n()) %>%
          arrange(desc(count))

`%nin%` <- Negate(`%in%`)


diff <- ko_count[ko_count$KO %nin% ko_map$KO,]
dim(diff)
write.table(diff,"data.enterprise/KAAS/ko_diff.ids", sep ="\t", quote = F,  row.names = F)
 
colnames(ko_ids)
colnames(ko_map)

kegg <- merge(ko_ids, ko_map, by = "KO", all.x = T ) 
dim(ko_ids)
dim(kegg)
colnames(kegg)

(kegg_category_count <- kegg %>%
              group_by(category) %>%
              summarise(count = n()) %>%
              arrange(desc(count)) %>%
              filter(category != "Human Diseases" & category != "-" ) %>%
              filter(!is.na(category)) %>%
              mutate(perc = round(100*count/sum(count),1)) %>%
              mutate(lab_pos = cumsum(perc) - 0.5* perc) %>%
              mutate(key = paste0(category, " [",perc,"%]")))

kegg_category_count$category <- factor(kegg_category_count$category, levels = rev(kegg_category_count$category))

kegg_doughnut  <-  ggplot(data = head(kegg_category_count, 10) , 
                  aes(x = 2,  y = perc,  fill = category)) +
                  geom_bar(width = 1,  stat = "identity") +
                  labs(x = NULL, y = NULL) +
                  guides(fill = guide_legend(reverse = TRUE)) +
                  geom_text_repel(aes(y = lab_pos,
                                      label = paste(perc,"%", sep = "")),
                                  min.segment.length = 0.75,
                                  col = "black",
                                  nudge_y = 0.05,
                                  nudge_x = 0.05,
                                  size = 6) +
                  theme(panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.title = element_text(size=18, colour = "black"),
                        axis.line = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(), 
                        legend.title = element_text(size=20, colour = "black"),
                        legend.text = element_text(size=20, colour = "black"),)+
                  coord_polar(theta = "y") +
                  xlim(0.5,2.5) +
                  scale_fill_simpsons(name = "Functional category")

KEGG_category <- kegg_category_count %>% select(category, count, perc)  
print(xtable(KEGG_category, digits=7, type = "latex"), file = "SuppInfoX/KEGG_sorted.tex", include.rownames = F)            

kegg_doughnut 
ggsave("FiguresX/KeggCategory.pdf", width = 16,  height = 8)

g <- ggplot_build(kegg_doughnut)
fill =  unique(g$data[[1]]["fill"])

fill <- rbind(fill, data.frame(fill = c("#000075","#f58231")))

kegg_top12 <- head(kegg_category_count, 10) %>%
              select(category)

kegg_top12 <- rbind(kegg_top12, 
      data.frame(category = c("Metabolism of other amino acids",
                              "Biosynthesis of other secondary metabolites")))

colour_map <- cbind(kegg_top12, fill)

colnames(kegg)
colnames(colour_map)

(kegg_name_count <- kegg %>%
                        group_by(KO ,name, category) %>%
                        summarise(count = n()) %>%
                        arrange(desc(count)) %>%
                        ungroup() %>%
                        mutate(label = paste0(name, " [",KO,"]")))

kegg_color_map  <- merge(kegg_name_count, colour_map, by = "category",  all.x = T) %>%
                   select(KO, name, category, label, count, fill) %>%
                   arrange(count) 

kegg_color_map$label  <- factor(kegg_color_map$label, levels = kegg_color_map$label)

tail(kegg_color_map)

kegg_barplot <- ggplot(data = kegg_color_map %>%  tail(30), aes( x = label,  y = count, fill= fill)) +
                   geom_col() +
                   coord_flip() +
                   scale_fill_identity() +
                   theme_bw() +
                   xlab("KEGG ortholog") +
                   ylab("Transcripts") +
                   theme( axis.text=element_text(size=16),
                        panel.grid = element_line(colour = "black",linetype="dashed",size=0.1),
                        axis.title = element_text(size=16),
                        legend.text = element_text(size=16)) +
                   scale_y_continuous(breaks = seq(0,60,10), limits = c(0,60))

kegg_barplot

ggsave("FiguresX/Keggnames.pdf", width = 12,  height = 8 )


###################################### PlantTFDB ################################

PlantTFDB <- read.delim("data.enterprise/PlantTFDB/TF_and_best1_in_Ath.list.txt", header = F)

colnames(PlantTFDB) <- c("TF_ID",	"Family","A_thaliana_BestHit", "Blast_evalue",	"Best_HitDescription")

PlantTFDB  <-  PlantTFDB %>%  
  filter(Blast_evalue <= 1e-05) %>%  
  arrange(desc(Blast_evalue))

colnames(PlantTFDB)

Family_TFDB <-  PlantTFDB %>% 
  group_by(Family) %>%
  summarise(Z_capensis =  n()) %>%
  arrange(desc(Z_capensis))


Zm_TF <- read.delim("data.enterprise/PlantTFDB/Zmr_TF_list.txt")

Family_Zm_TF <- Zm_TF %>% 
  group_by(Family) %>%
  summarise(Z_marina =  n()) %>%
  arrange(desc(Z_marina))

Family_TF <- merge(Family_Zm_TF, Family_TFDB, all.x = T, all.y = T) 

Family_TF[is.na(Family_TF)] <- 0

tmp <- list()
tmp$sums <- rowSums(Family_TF[c("Z_marina","Z_capensis")])

Family_TF <- Family_TF[order(tmp$sums, decreasing = F),]

Family_TF_melt <- melt(Family_TF)

Family_TF_melt$Family <- factor(Family_TF_melt$Family,  levels = Family_TF$Family)

FamilyTF <- ggplot(Family_TF_melt, aes(x = Family, y = value, fill = variable)) + 
                geom_bar(stat = "identity", position = position_dodge()) +
                coord_flip()+
                labs(x = "TF Family", y ="Transcription factors" ) +
                theme( panel.grid = element_line(colour = "black",linetype="dashed",size=0.1),
                      panel.background = element_blank(),
                      axis.text = element_text(size=12, colour = "black"),
                      axis.title = element_text(size=14, colour = "black"),
                      legend.position = "top",
                      legend.text = element_text(size=12, colour = "black"),
                      legend.background = element_rect(colour = "black", fill=NA, size=0.5),
                      strip.background = element_rect(colour = "black", fill=NA, size=1),
                      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
                scale_fill_manual(name = "Species", values =c("#808000","#00B7B1"), labels = c(expression(italic("Zostera marina")), expression(italic("Zostera capensis")))) +
                ylim(c(0,200))

FamilyTF

ggsave("FiguresX/Family_TF_melt.pdf")  

print(xtable(arrange(Family_TF, desc(Z_capensis)), digits=7, type = "latex"), file = "SuppInfoX/TF.tex", include.rownames = F) 

############################# Swiss-Prot GO ####################################

Regulation_ref <- read.table("data.enterprise/Blast_top_hits/Regulation_ref.tsv", header = T)


(Regulation <- Regulation_ref %>%
  group_by(DEG) %>%
  summarise(count = n()))

DEG_count  <-  ggplot(Regulation, aes(x = DEG, y= count, fill = DEG)) +
                    geom_col() +
                    coord_flip() +
                    ylab("Number of genes") +
                    xlab("Response") +
                    theme( panel.grid = element_line(colour = "black",linetype="dashed",size=0.05),
                         panel.background = element_blank(),
                         axis.text = element_text(size=16, colour = "black"),
                         axis.title = element_text(size=16, colour = "black"),
                         legend.text = element_text(size=16, colour = "black"),
                         legend.background = element_rect(colour = "black", fill=NA, size=0.5),
                         strip.background = element_rect(colour = "black", fill=NA, size=1),
                         panel.border = element_rect(colour = "black", fill=NA, size=1), 
                         legend.position = 'none') +
                         scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
                         scale_fill_manual(values=c("#1F77B4FF","#FF7F0EFF"))
DEG_count

ggsave("FiguresX/RegulationX.pdf", height = 40 , width = 210, units = "mm")


(uniprot_ref <- read.delim("data.enterprise/Blast_top_hits/uniprot_sprot.ref", sep="\t", header = F) %>% 
    `colnames<-`(c("ACC","name"))) %>%
head()

colnames(Regulation_ref)
colnames(uniprot_ref)
dim(Regulation_ref)
(top_DEG_names <- Regulation_ref %>% pull(ACC))

DEG_data <- merge(Regulation_ref, uniprot_ref, by = "ACC", all.x = T) %>%
            dplyr::select(ACC, name, log2FoldChange, lfcSE, padj) %>%
            arrange(desc(log2FoldChange))


gene_names <- read.csv("data.enterprise/Blast_top_hits/DEG_genames.tab", sep = "\t", header = T)
colnames(gene_names)[3] <- "Gene_name"

gene_symbols <- gene_names %>%
            dplyr::select(Entry, Gene_name) %>%
            rename(ACC = Entry)

colnames(gene_symbols)

DEG_datax <- merge(DEG_data, gene_symbols) %>%
             dplyr::select(ACC, name, Gene_name, log2FoldChange, lfcSE, padj)

colnames(DEG_datax)

dim(DEG_datax)

head(DEG_datax)

print(xtable(DEG_datax, digits= 12, type = "latex"), file = "SuppInfoX/DEG_data.tex", include.rownames = F)

################################## Volcano plot #########################################
#show_col(pal_d3("category10")(10))

(dseq2_results <- read.csv("data.enterprise/Blast_top_hits/TrX_padj1.tsv", sep = "\t") %>%
   drop_na() %>%
   mutate(DEG_type = case_when(padj <= 0.05 &  log2FoldChange >= 1 ~ 'Upregulated',
                               padj <= 0.05 &  log2FoldChange < -1 ~ 'Downregulated',
                               padj > 0.05  |   log2FoldChange < 1 ~ 'Non-significant')))

gene_names <- read.csv("data.enterprise/Blast_top_hits/DEG_genames.tab", sep = "\t", header = T)
colnames(gene_names)[3] <- "gene_names"

dim(dseq2_results)
dseq2_volcono <-  merge(dseq2_results, gene_names, by.x="SwissProt", by.y = "Entry",  all.x = T)
dim(dseq2_results)

h_lfc <- dseq2_volcono %>%  drop_na() %>% arrange(log2FoldChange) %>% head()
t_lfc <- dseq2_volcono %>%  drop_na() %>% arrange(log2FoldChange) %>% tail()
t_padj <- dseq2_volcono %>%  drop_na() %>% filter(padj < 1e-3) 

dseq2_results_label <-rbind(h_lfc, t_lfc, t_padj) %>%  distinct()

dseq2_results_label[dseq2_results_label$gene_names =="",]$gene_names <- paste0(dseq2_results_label[dseq2_results_label$gene_names =="",]$Entry.name,"*")

dseq2_results_label %>%
  dplyr::select(SwissProt, log2FoldChange, padj,DEG_type, gene_names, Organism)

summary(dseq2_results$log2FoldChange)
summary(dseq2_results$padj)

dseq2_volcono$DEG_type <- factor(dseq2_volcono$DEG_type, levels = c("Upregulated","Downregulated","Non-significant"))

#show_col(c("#1F77B4FF", "#808080","#FF7F0EFF"))

volcano <-ggplot(dseq2_volcono , aes(y = -log10(padj), x = log2FoldChange, alpha =.9)) + 
  geom_point(size = 3, aes(color = DEG_type) ) + 
  scale_color_manual(values=c("#FF7F0EFF","#1F77B4FF","#808080"), name = 'Gene expression') +
  ggrepel::geom_text_repel(data = dseq2_results_label, 
                           aes(y = -log10(padj), 
                               x = log2FoldChange, 
                               label = gene_names), 
                           size = 3.5) +
  xlab(expression(paste(Log[2]," ", fold," ", change))) +
  ylab(expression(paste(-Log[10],"(P-value)"))) +
  theme_bw() +
  geom_hline(yintercept= -log10(0.05), size = 0.1, colour="#00FF00", linetype="dashed") +
  geom_vline(xintercept=c(1,-1), size = 0.1, colour=c("#BB0000","#1F77B4FF"), linetype="dashed") +
  theme(legend.background = element_rect(colour = "darkgray"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16, vjust = 1),
        legend.title=element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "top",
        panel.grid = element_blank())  +
  scale_alpha(guide = 'none') +
  scale_x_continuous(limits = c(-13,7), breaks = seq(-14,8,2)) +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2.5))

volcano      

ggsave("FiguresX/volcano.pdf", width = 210, height = 180, units = c("mm"))


############################## Regulation GO term ##############################
colnames(Swissprot_GO_slim)
dim(Swissprot_GO_slim)
GO_ref <-  Swissprot_GO_slim %>%
           dplyr::select(-transcript) %>%
           distinct() %>%
           dplyr::rename(ACC=Uniprot_ID )
           
length(unique(GO_ref$GO_accession))

dim(GO_ref)
colnames(GO_ref)

(Regulation_ref <- read.table("data.enterprise/Blast_top_hits/Regulation_ref.tsv", header = T) %>%
                   dplyr::select(ACC, DEG)) %>%
                   head(10)

dim(Regulation_ref)

GO_regulation <-merge(Regulation_ref, GO_ref, by = "ACC", all.x = T)
dim(GO_regulation)

GO_regulation %>%
     select(ACC, DEG, GO_accession, goslim_domain, term) %>%
     dplyr::filter(DEG == "Downregulated") %>%
     dplyr::filter(goslim_domain == "biological_process") %>%
     dplyr::filter(!is.na(GO_accession))

noquote(colnames(GO_regulation))

DEGGO_data <-  GO_regulation %>% 
               dplyr::select(ACC,DEG,GO_accession,goslim_domain,goslim_term)

print(xtable(DEGGO_data, digits=7, type = "latex"), file = "SuppInfoX/DEGGO_data.tex", include.rownames = F)            

(DEG_no_GO <- GO_regulation[is.na(GO_regulation$GO_accession),] %>%
              dplyr::select(ACC, DEG) %>%
              arrange(DEG))

my_deg <- DEGGO_data[!is.na(DEGGO_data$GO_accession),]


my_deg %>% 
  select(ACC, DEG, GO_accession, term, goslim_domain) %>% 
  filter(DEG == "Upregulated") %>%
  filter(DEG == "P42730")


unique(my_deg$ACC) %T>%
           print() %>%
           length()

DEG_no_GO %>%
         group_by(DEG) %>%
         summarise(count =  n())


GO_regulation_count <- GO_regulation   %>%  
                       filter(!is.na(GO_accession)) %>%
                       group_by(GO_accession, DEG, goslim_domain, goslim_term) %>%
                       summarise(count = n()) %>%
                       arrange(desc(count)) 


GO_regulation_count %>%
                    filter(goslim_domain == "molecular_function")


head(GO_regulation_count)


GO_regulation_count <- GO_regulation_count %>%
                     arrange(goslim_domain, desc(count))

unique(GO_regulation_count$goslim_term)
     
GO_regulation_count$goslim_term <- factor(GO_regulation_count$goslim_term, levels =  unique(GO_regulation_count$goslim_term))

GO_regulation_count$goslim_domain <- factor(GO_regulation_count$goslim_domain,  
                                            levels = c("biological_process", "molecular_function","cellular_component"),
                                            labels = c("Biological process", "Molecular function","Cellular component"))

write.table(GO_regulation_count, file = "data.enterprise/Blast_top_hits/GO_regulation_count.tsv", row.names = FALSE, quote = FALSE, sep ="\t", col.names = T)

GO_reg <- ggplot(GO_regulation_count, aes(x = goslim_term ,y = count, fill = DEG)) +
               geom_col(position = position_dodge2(width = 1, preserve = "single")) + 
               coord_flip() +
               ylab("Number of terms") +
               facet_grid(rows = vars(goslim_domain), scales = "free_y", switch = "y", space = "free_y", labeller = label_wrap_gen(10)) +
            theme(panel.grid = element_line(colour = "black",linetype="dashed",size=0.05),
                  panel.background = element_blank(),
                  axis.text = element_text(size=14, colour = "black"),
                  axis.title = element_text(size=14, colour = "black"),
                  legend.position = "top",
                  legend.title = element_text(size=14, colour = "black"),
                  legend.text = element_text(size=14, colour = "black"),
                  legend.background = element_rect(colour = "black", fill=NA, size=0.1),
                  strip.background = element_rect(colour = "black", fill=NA, size=1),
                  panel.border = element_rect(colour = "black", fill=NA, size=1),
                  strip.placement = "outside",
                  axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
                  axis.title.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  strip.text = element_text(size=14, colour = "black"), 
                  strip.background.x=element_rect(color = NA,  fill=NA), 
                  strip.background.y=element_rect(color = "black",  fill=NA)) + 
                  scale_y_continuous(breaks = seq(0, 14, 2), limits = c(0,14))+
                  scale_fill_d3(name = "Gene response")

GO_reg


ggsave("FiguresX/DEG_GO_count.pdf", width = 225, height  = 200, units = c("mm"))

GO_regulation_count %>% data.frame()



####################################### External Data #################################

(stress_proteins <- read.csv("data.enterprise/ExternalX/external_study2.tsv",  sep = "\t",  na.strings = c("","NA")))
colnames(stress_proteins)[1] <- "Accession"
colnames(stress_proteins)
StressProteins_blastx   <- read.csv("data.enterprise/ExternalX/AllProteins.blastx", sep = "\t") %>%
                         `colnames<-`(c("qseqid1", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend","sstart","send", "evalue" ,"bitscore"))
colnames(stress_proteins)
colnames(StressProteins_blastx)

head(stress_proteins)
head(StressProteins_blastx)

lapply(list(stress_proteins, StressProteins_blastx), dim)
lapply(list(stress_proteins, StressProteins_blastx), colnames)

(stress_proteins_results <- merge(stress_proteins, StressProteins_blastx, by.x = "Entry",  by.y = "qseqid1", all.x = TRUE))

`%nin%` <- Negate(`%in%`)

stress_proteins[stress_proteins$Entry %nin% StressProteins_blastx$qseqid1,] 

(GOI_data <- stress_proteins_results %>%
                       select(Accession, Gene.names, Protein.names, sseqid, pident, length,  evalue) %>%
                       arrange(Accession))

print(xtable(GOI_data, digits=7, type = "latex"), file = "SuppInfoX/GOI_data.tex", include.rownames = F)            

################################## DGE Gene-set ################################

DGE_genset <- read.csv("data.enterprise/DGE/gene_set.txt",  sep = "\t", header = F) %>%
              `colnames<-`(c("Uniprot_ID"))
dim(DGE_genset)

Blast_Swissprot  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.ids", sep = "\t") %>%
                   `colnames<-`(c("transcript","Uniprot_ID"))

DGE_SwissProt <- merge(DGE_genset, Blast_Swissprot, by.x = "Uniprot_ID", all.x = T)
length(DGE_genset$Uniprot_ID) == length(unique(DGE_SwissProt$Uniprot_ID))

DGE_SwissProt_uniq <- DGE_SwissProt %>%
                      group_by(transcript, Uniprot_ID) %>%
                      summarise(count = n()) %>%
                      arrange(desc(count))

DGE_SwissProt_uniq %>% 
                      group_by(Uniprot_ID) %>%
                      summarise(count = n()) %>%
                      dim()

write.table(DGE_SwissProt_uniq, 
                      file = "data.enterprise/DGE/DGE_SwissProt_uniq.tsv",
                      row.names = FALSE, 
                      quote = FALSE, 
                      sep ="\t", 
                      col.names = T)

(DGE_transcript_uniq <- DGE_SwissProt_uniq %>% 
                        group_by(transcript) %>%
                        summarise(count = n())) 

Swissprot_GO_map  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.map", sep = "\t") %>%
                     `colnames<-`(c("transcript","Uniprot_ID","GO_accession","domain","term"))
dim(Swissprot_GO_map)

lapply(list(DGE_SwissProt, Swissprot_GO_map), colnames)
lapply(list(DGE_SwissProt, Swissprot_GO_map), dim)
DGE_GO_SwissProt <- merge(DGE_genset, Swissprot_GO_map, by = "Uniprot_ID", all.x = T)
dim(DGE_GO_SwissProt)
colnames(DGE_GO_SwissProt)
 
`%nin%` <- Negate(`%in%`)
DGE_genset[DGE_genset$Uniprot_ID %nin% Swissprot_GO_parsed$Uniprot_ID,]
dim(DGE_GO_SwissProt)

plant_GO_slim <- read.csv("data.enterprise/Blast_top_hits/goslim_plant.tsv", sep = "\t",  header = F) %>%
                `colnames<-`(c("GO_accession", "goslim_domain", "goslim_term"))
 colnames(DGE_GO_SwissProt)

DGE_plant <- DGE_GO_SwissProt[DGE_GO_SwissProt$GO_accession %in% plant_GO_slim$GO_accession,]
DGE_GOPlant <- merge(DGE_GO_SwissProt,plant_GO_slim)
lapply(list(DGE_GOPlant, DGE_plant), dim)
colnames(DGE_GOPlant)

DGE_GOPlant %>%
   select(Uniprot_ID,transcript,domain)



DGE_count <- DGE_GOPlant %>%
             group_by(goslim_domain) %>%
             summarise(count = n())
            

################################## Expasy-Enzyme ###############################
Blast_Swissprot  <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.ids", sep = "\t") %>%
                    `colnames<-`(c("transcript","Uniprot_ID"))
dim(Blast_Swissprot)



enzyme_uniprot_map <- read.csv("data.enterprise/Blast_top_hits/enzyme_uniprot.map",sep = "\t", header = F) %>%
                      `colnames<-`(c("Uniprot_ID","EC_number"))
dim(enzyme_uniprot_map)

Blast_enzyme <- merge(Blast_Swissprot,  enzyme_uniprot_map, by = "Uniprot_ID")
colnames(Blast_enzyme)
dim(Blast_enzyme) 


################################### HSP ########################################

HSP_name <- read.csv("data.enterprise/Blast_top_hits/HSP.names", sep = "\t") %>%
            `colnames<-`(c("UniprotID","name"))


Transcript2Uniprot <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.ids", sep = "\t") %>%
            `colnames<-`(c("Transcript","UniprotID"))

Blast_results <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.blastx.uniprot", sep = "\t") %>%
                `colnames<-`(c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart" ,"q_end","sstart","send","evalue","bitscore"))

HSP_transcripts <- merge(Transcript2Uniprot, HSP_name) %>%
             select(Transcript, UniprotID, name)

HSP_ref  <- HSP_transcripts %>%
            select(UniprotID, name) %>%
            distinct()


HSP_Blast_results <- Blast_results[Blast_results$sseqid  %in% HSP_transcripts$UniprotID,] 

`%nin%` <- Negate(`%in%`)

length(HSP_transcripts$Transcript) ==  length(HSP_Blast_results$qseqid)
lapply(list(HSP_transcripts, HSP_Blast_results), dim)

hsp_data <- merge(HSP_Blast_results,  HSP_ref,  by.x = "sseqid", by.y = "UniprotID") %>%
            group_by(sseqid, name) %>%
            summarise(count = n()) %>%
            arrange(desc(count))


write.table(hsp_data, 
            file = "data.enterprise/Blast_top_hits//hsp_data.tsv",
            row.names = FALSE, 
            quote = FALSE, 
            sep ="\t", 
            col.names = T)

print(xtable(hsp_data, digits=7, type = "latex"), file = "SuppInfoX/HSP_data.tex", include.rownames = F)

################################################################################

