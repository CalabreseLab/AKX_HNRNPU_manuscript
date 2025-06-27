############## Silhouette Coefficient
################################
library(igraph)
library(graphkernels)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(cluster)


dt<-fread('comb_allgenes_cor_08092024_long.csv',header=T)

xist<-'Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3'
airn<-'Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced'
ot1<-'Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced'

allcom<-fread('ldcommunity_allgenes_08092024.csv',header=T)
# 13832, actually 13831 genes, as there are one ID columns
# 5464 genes that contains NA in the Weights col

silo<-function(gene1,dt,allcom) {
  
  # Filter rows for each gene using data.table's syntax
  dt1 <- dt[gene == gene1]
  dt1[, gene := NULL]
  setnames(dt1, old = "Weight", new = "weight")
  
  # filter and keep only the positive edges
  dt1_pos<-dt1[weight>=0]
  
  # convert the similarity to distance for this calculation
  dt1_pos$weight<-(1-dt1_pos$weight)
  
  # Create the igraph object directly from the data.table
  # as the column name is weight it will be recognized as edge weight in the graph
  graph1 <- graph_from_data_frame(d = dt1_pos, directed = FALSE)
  # check edges E(graph1)
  # check weights
  # E(graph1)$weight
  # check nodes V(graph1)
  
  # Compute the shortest path distance matrix
  dist1 <- distances(graph1,weights=E(graph1)$weight)
  
  # check if there are Inf in the distance, happens when there are isolated nodes
  # that are not connected to the main graph
  if (any(is.infinite(dist1))) {
    # return NA if there are isolated nodes
    # can also filter out the isolated nodes 
    # Extract the Largest Connected Component
    return(NA)
    
  } else {
    # get community assignment
    com1<-allcom[[gene1]]
    names(com1)<-allcom$Id
    
    # Reorder community labels to match the distance matrix order
    com1 <- com1[rownames(dist1)]
    
    # if each component is by itself a community
    if (max(com1)==length(com1)) {
      return(NA)
    } else {
      
      # Compute the Silhouette Coefficient
      sil <- silhouette(com1, dist1)
      
      # Get the average Silhouette width
      avg_sil_width <- mean(sil[, 'sil_width'])
      
      return(avg_sil_width)
      
    }
    
  }
  
}


xist_s<-silo(xist,dt,allcom)
airn_s<-silo(airn,dt,allcom)
ot1_s<-silo(ot1,dt,allcom)

allgene<-colnames(allcom)[-1]
allgene_s<-vector(mode='numeric',length=length(allgene))

for (g in 1:length(allgene)) {
  
  if(g %% 1000 == 0){print(g)}
  
  gene1<-allgene[g]
  
  allgene_s[g]<-silo(gene1,dt,allcom)
  
}

allgene_s<-as.data.frame(allgene_s)
colnames(allgene_s)<-'silo'
allgene_s$gene<-allgene

xquant<-round(sum(xist_s>allgene_s$silo,na.rm=T)/sum(!is.na(allgene_s$silo)),digits=4)*100
aquant<-round(sum(airn_s>allgene_s$silo,na.rm=T)/sum(!is.na(allgene_s$silo)),digits=4)*100
kquant<-round(sum(ot1_s>allgene_s$silo,na.rm=T)/sum(!is.na(allgene_s$silo)),digits=4)*100


xak<-data.frame(gene=c(paste0('Xist(',xquant,'%)'),
                       paste0('Airn(',aquant,'%)'),
                       paste0('Kcnq1ot1(',kquant,'%)')),
                silo=c(xist_s,airn_s,ot1_s))


hbreaks<-seq(from=-0.05, to=0.75, by=0.025)

sum(!is.na(allgene_s$silo))

p<-ggplot() +
  geom_histogram(data = allgene_s, aes(x = silo, y = after_stat(count)),  fill='lightgrey',color='darkgrey', breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=2000, by=500)) +
  labs(x = "Mean Silhouette Width") +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24))+ 
  geom_vline(data = xak, aes(xintercept = silo), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = xak, aes(x = silo, y = c(1000,1250,1500), label = gene), 
                  size = 8, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")



ggsave("RBP_XAK_hist_Silhouette.pdf", plot = p, width = 5, height = 4, units = "in")

