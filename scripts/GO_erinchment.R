library(ViSEAGO)
library(dplyr)
library(topGO)


`%nin%` <- Negate(`%in%`)
Background_SwissProt <- read.csv("data.enterprise/Blast_top_hits/EvigeneX.swissprot", header = T) %>% `colnames<-`(c("ACC"))
dim(Background_SwissProt)
DGE_SwissProt <- read.csv("data.enterprise/Blast_top_hits/DGE.ids", header = F) %>% `colnames<-`(c("ACC"))
dim(DGE_SwissProt)
############################### GO terms BP ####################################
GO_plant <- read.csv("data.enterprise/Blast_top_hits/Swissprot_GO_slim.tsv", sep = "\t") %>%
                      filter(goslim_domain == "biological_process") %>%
                      dplyr::select(transcript,Uniprot_ID,GO_accession)

GO_plant$taxid <- "Z.capensis"
GO_plant$evidence <- NA
colnames(GO_plant)
GO_plant <- GO_plant %>%
            `colnames<-`(c("gene_symbol","gene_id","GOID","taxid","evidence"))
colnames(GO_plant)

############################### Custom2GO ######################################

Custom <- ViSEAGO::Custom2GO("data.enterprise/Blast_top_hits/GO_plant.tsv")

ViSEAGO::available_organisms(Custom)

myGENE2GO <-ViSEAGO::annotate(
  "Z.capensis",
  Custom
)

####################################### TopGO ###################################
p_value = 0.05
nodeSize = 5

BP_viSEAGO  <-ViSEAGO::create_topGOdata(
  geneSel=DGE_SwissProt$ACC,
  allGenes=Background_SwissProt$ACC,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=nodeSize
)


# perform topGO tests
BP_topGO <- topGO::runTest(
  BP_viSEAGO,
  algorithm ="classic",
  statistic = "fisher",
  cutOff=p_value
)

BP_Results <- ViSEAGO::merge_enrich_terms(
  cutoff=p_value,
  Input=list(
    condition=c("BP_viSEAGO","BP_topGO")
  )
)

ViSEAGO::show_table(BP_Results)

myBP_SS<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_Results
)

# compute all available Semantic Similarity (SS) measures
myGO_BP <-ViSEAGO::compute_SS_distances(
  myBP_SS,
  distance="Wang"
)

ViSEAGO::MDSplot(myGO_BP)


Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGO_BP,
  showIC=FALSE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)###################           
# Display the clusters-heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

###################    
# print the clusters-heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  type = "GOterms",
  file = "SuppInfoX/cluster_heatmap_Wang_wardD2.pdf",
  plotly_update = T
)

#SwissProt2Entrez <- read.csv("data.enterprise/Blast_top_hits/Evigene.entrez", sep = "\t", header = T) 
#dim(SwissProt2Entrez)
#NotMappedSwissprot <- Background_SwissProt[Background_SwissProt$ACC  %nin%  SwissProt2Entrez$From,]
#length(SwissProt2Entrez$From)
#SwissProt2Entrez[duplicated(SwissProt2Entrez$From),]
# length(unique(SwissProt2Entrez$From))
# length(NotMappedSwissprot)
# length(unique(NotMappedSwissprot))
# ################################################################################

# DGE_SwissProt2Entrez <- read.csv("data.enterprise/Blast_top_hits/DGE.entrez", sep = "\t", header = T) 
# dim(DGE_SwissProt2Entrez)
# NotMappedDGE <- DGE_SwissProt[DGE_SwissProt$ACC  %nin%  DGE_SwissProt2Entrez$From,] 
# length(NotMappedDGE)
# 
# Bioconductor<-ViSEAGO::Bioconductor2GO()
# 
# myGENE2GO <-ViSEAGO::annotate(
#   "org.At.tair.db",
#   Bioconductor
# )
# 
# myGENE2GO

