# for each gene calculate whether the all double/triple is intra or inter community
# add the status of inter/intra community of XAK to the double/triple table
# for the intra community protein double/triple pair check the adjpval and see if prev/rare
# take the ratio of uncommon/common pairs for each gene
# check where XAK ratio stand among all genes


allgenes<-read.csv('genelist13831.csv',header=T)
allgenes<-allgenes$x

allfiles<-dir('/work/users/s/h/shuang9/rip/multicov_ld_nodes/')
ldnames<-gsub('_nodes_ld.csv','',allfiles)
allfiles<-allfiles[ldnames %in% allgenes]

dbl<-read.csv('alldoublepro_ld_pval_adj_06182025_raw.csv',header=T)

# initiate dataframe to store the whether each protein pair is inter-0 or intra-1


dbl_commu<-matrix(data=NA,nrow=length(allgenes),ncol=nrow(dbl),)
rownames(dbl_commu)<-allgenes
colnames(dbl_commu)<-dbl$proset

dbl_pro<-strsplit(dbl$proset,'_',fixed=T)
dbl_pro <- do.call(rbind, dbl_pro)

# loop thru genes

for (n in 1:length(allgenes)) {
  
  if(n %% 100==0) {print(n)}
  
  tf<-allfiles[n]
  
  tnodes<-read.csv(paste0('../multicov_ld_nodes/',tf),header=T)
  
  # go thru double pairs
  for (m in 1:nrow(dbl_pro)) {
    
    p1<-tnodes$Color[which(tnodes$Id==dbl_pro[m,1])]
    p2<-tnodes$Color[which(tnodes$Id==dbl_pro[m,2])]
    # same community-1 or dif commmu-0
    if (p1==p2) {
      dbl_commu[n,m]<-1
    } else {
      dbl_commu[n,m]<-0
    }
    
  }
  

}

write.csv(dbl_commu,'allgenes_doublepro_community.csv')

# add XAK to the tables
dbl$Xist<-as.numeric(dbl_commu['Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3',])

dbl$Airn<-as.numeric(dbl_commu['Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced',])

dbl$Kcnq1ot1<-as.numeric(dbl_commu['Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced',])


write.csv(dbl,'alldoublepro_ld_pval_adj_06182025_raw.csv',row.names = F)

###############################################################
# calculate the prev and rare ratio for all protein pairs that are intra-1


dbl<-read.csv('alldoublepro_ld_pval_adj_06182025_raw.csv',header=T)
dbl_commu<-read.csv('allgenes_doublepro_community.csv',header=T,row.names = 1)

dbl_prev<-which(dbl$prev_adjp<0.05)

dbl_rare<-which(dbl$rare_adjp<0.05)


allgenes<-row.names(dbl_commu)

dbl_stats<-data.frame(genes=allgenes,
                      total_intra=0,
                      rare=0,
                      prev=0)


dbl_stats$total_intra<-rowSums(dbl_commu)


# rare edges are the ones with sig rare_adjp
# sig p value edges are the ones that appears in more genes than expected
dbl_stats$rare<-rowSums(dbl_commu[,dbl_rare])

dbl_stats$prev<-rowSums(dbl_commu[,dbl_prev])

dbl_stats$rare_ratio<-dbl_stats$rare/dbl_stats$total_intra

write.csv(dbl_stats,'allgenes_double_intra_rare_ratio_06192025.csv',row.names = F)


#there are NAs in each stats -- no total intra edges
sum(is.na(dbl_stats$rare_ratio)) # 1
#sum(is.na(tpl_stats$uncommon_ratio)) # 1
# 7230
# Mal2_chr15_54571192_54602847_+_ENSMUSG00000024479.3_ENSMUSG00000024479.3.unspliced
# this one has each node into one community, 27 communities in total so no intra edges

# remove the NA row to be excluded from the next analysis
dbl_stats<-dbl_stats[!is.na(dbl_stats$rare_ratio),]


xidx<-which(dbl_stats$genes=='Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3')
aidx<-which(dbl_stats$genes=='Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced')
kidx<-which(dbl_stats$genes=='Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced')

dbl_xak<-dbl_stats[c(xidx,aidx,kidx),]
#tpl_xak<-tpl_stats[c(xidx,aidx,kidx),]

dbl_xak$genename<-c('Xist','Airn','Kcnq1ot1')

dbl_xak$total_intra_pct<-0
dbl_xak$rare_pct<-0
dbl_xak$rare_ratio_pct<-0

# add percentile for rare, total intra and rare ratio
for (n in 1:nrow(dbl_xak)){
  
  dbl_xak$total_intra_pct[n]<-sum(dbl_xak$total_intra[n]>dbl_stats$total_intra)*100/nrow(dbl_stats)
  
  dbl_xak$rare_pct[n]<-sum(dbl_xak$rare[n]>dbl_stats$rare)*100/nrow(dbl_stats)
  
  dbl_xak$rare_ratio_pct[n]<-sum(dbl_xak$rare_ratio[n]>dbl_stats$rare_ratio)*100/nrow(dbl_stats)
  
}

dbl_xak$total_intra_label<-paste0(dbl_xak$genename,'(',round(dbl_xak$total_intra_pct,2),'%)')
dbl_xak$rare_label<-paste0(dbl_xak$genename,'(',round(dbl_xak$rare_pct,2),'%)')
dbl_xak$rare_ratio_label<-paste0(dbl_xak$genename,'(',round(dbl_xak$rare_ratio_pct,2),'%)')


library(ggplot2)
library(ggrepel)

hbreaks<-seq(from=0, to=1, by=0.05)

p<-ggplot() +
  geom_histogram(data = dbl_stats, aes(x = rare_ratio, y = after_stat(count)),  fill='lightgrey',color='darkgrey', breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=2400, by=500)) +
  labs(x = "Rare/Total Intracommunity Edges",y = "Count") +
  coord_cartesian(xlim=c(0,0.8))+
  scale_x_continuous(breaks=seq(from=0, to=0.8, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  theme(panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24))+ 
  geom_vline(data = dbl_xak, aes(xintercept = rare_ratio), color = "#1b9e77", linetype = "solid",linewidth=1)+
  geom_text_repel(data = dbl_xak, aes(x = rare_ratio, y = c(1500,1800,2200), label = rare_ratio_label), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("allgenes_double_rare_ratio_XAK_hist_06192025.pdf", plot = p, width = 7.5, height = 5, units = "in")


hbreaks<-seq(from=0, to=100, by=5)

p<-ggplot() +
  geom_histogram(data = dbl_stats, aes(x = rare, y = after_stat(count)),  fill='lightgrey',color='darkgrey', breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=3000, by=500)) +
  labs(x = "Rare Edges",y = "Count") +
  coord_cartesian(xlim=c(0,80),ylim=c(0,3000))+
  scale_x_continuous(breaks=seq(from=0, to=80, by=20))+
  theme(panel.background=element_rect(fill='white'),
        #plot.margin = margin(2, 2, 2, 2, "pt"),
        # panel.grid.major=element_line(color='grey',linewidth =0.3),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24))+ 
  geom_vline(data = dbl_xak, aes(xintercept = rare), color = "#1b9e77", linetype = "solid",linewidth=1)+
  geom_text_repel(data = dbl_xak, aes(x = rare, y = c(1800,2300,2800), label = rare_label), 
                  size = 8, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)

ggsave("allgenes_double_rare_XAK_hist_06192025.pdf", plot = p, width = 7.5, height = 5, units = "in")



