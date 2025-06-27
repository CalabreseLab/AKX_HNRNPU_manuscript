#################################################################
# calculate modularity of XAK network
library(igraph)
library(data.table)
library(ggplot2)
library(ggrepel)

xist<-'Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3'
airn<-'Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced'
ot1<-'Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced'

# this is the community assignment
allcom<-fread('ldcommunity_allgenes_08092024.csv',header=T)
# 13832, actually 13831 genes, as there are one ID columns
# 5464 genes that contains NA in the Weights col

dt<-fread('comb_allgenes_cor_08092024_long.csv',header=T)
# 6772545/19295=351

modu<-function(gene1,dt,allcom) {
  
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
  
  # get community assignment
  com1<-allcom[[gene1]]
  names(com1)<-allcom$Id
  
  # Reorder community labels to match the V(g) order
  com1 <- com1[V(graph1)$name]
  
  # Compute modularity
  Q <- modularity(graph1, com1)
  
  return(Q)
  
  
}


xist_mod<-modu(xist,dt,allcom) # 0.1746962
airn_mod<-modu(airn,dt,allcom) # 0.06132606
ot1_mod<-modu(ot1,dt,allcom) # -0.02468401

allgene<-colnames(allcom)[-1]
allgene_mod<-vector(mode='numeric',length=length(allgene))

for (g in 1:length(allgene)) {
  
  if(g %% 1000 == 0){print(g)}
  
  gene1<-allgene[g]
  
  allgene_mod[g]<-modu(gene1,dt,allcom)
  
}

allgene_mod<-as.data.frame(allgene_mod)
colnames(allgene_mod)<-'modularity'
allgene_mod$gene<-allgene

xquant<-round(sum(xist_mod>allgene_mod$modularity,na.rm=T)/sum(!is.na(allgene_mod$modularity)),digits=4)*100
# 98.89
aquant<-round(sum(airn_mod>allgene_mod$modularity,na.rm=T)/sum(!is.na(allgene_mod$modularity)),digits=4)*100
# 82.96
kquant<-round(sum(ot1_mod>allgene_mod$modularity,na.rm=T)/sum(!is.na(allgene_mod$modularity)),digits=4)*100
# 11.19

xak<-data.frame(gene=c(paste0('Xist(',xquant,'%)'),
                       paste0('Airn(',aquant,'%)'),
                       paste0('Kcnq1ot1(',kquant,'%)')),
                modularity=c(xist_mod,airn_mod,ot1_mod))

range(allgene_mod$modularity)

hbreaks<-seq(from=-0.05, to=0.3, by=0.01)

sum(!is.na(allgene_mod$modularity))


p<-ggplot() +
  geom_histogram(data = allgene_mod, aes(x = modularity, y = after_stat(count)),  fill='lightgrey',color='darkgrey', breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=1500, by=500)) +
  labs(x = "Modularity") +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24))+ 
  geom_vline(data = xak, aes(xintercept = modularity), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = xak, aes(x = modularity, y = c(1000,1250,1500), label = gene), 
                  size = 8, nudge_x = 0.01, direction = "y", hjust = -0.1, vjust = -0.5,
                  segment.color = "darkgray",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid")



ggsave("RBP_XAK_hist_Modularity.pdf", plot = p, width = 5, height = 4, units = "in")

