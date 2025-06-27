#Want to do biotype analysis first on the top 1000 transcripts similar to X, K, A by RBP binding intensity

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

comb<-read.csv('allgenes_XAKRIPr_XAKnetworkr_common_sil_modu_biotype_06242025.csv',header=T)

comb<-comb[,c(1:4,14)]

# edit biotype
comb$biotype[which(comb$biotype!='protein_coding' & comb$biotype!='lncRNA')]<-'others'

# get top1000 xak, exclude XAK from the list
top1000_xist <- comb[order(comb$Xist_RIP_r,decreasing = TRUE)[2:1001],]
top1000_airn <- comb[order(comb$Airn_RIP_r,decreasing = TRUE)[2:1001],]
top1000_ot1 <- comb[order(comb$Kcnq1ot1_RIP_r,decreasing = TRUE)[2:1001],]
  

biotype_counts_xist <- as.data.frame(table(top1000_xist$biotype))
biotype_counts_airn <- as.data.frame(table(top1000_airn$biotype))
biotype_counts_ot1 <- as.data.frame(table(top1000_ot1$biotype))

biotype_counts_chrom <- as.data.frame(table(comb$biotype))

# Rename the columns
colnames(biotype_counts_xist) <- c("Biotype", "Count_xist")
colnames(biotype_counts_ot1) <- c("Biotype", "Count_ot1")
colnames(biotype_counts_airn) <- c("Biotype", "Count_airn")
colnames(biotype_counts_chrom) <- c("Biotype", "Count_chrom")

joined_data <- full_join(biotype_counts_chrom, biotype_counts_xist, by = "Biotype") %>%
  full_join(biotype_counts_airn, by = "Biotype")%>%
  full_join(biotype_counts_ot1, by = "Biotype")


#use the code below to generate spliced and unspliced info for each biotype
dt_x<-as.data.table(top1000_xist)[,.N,by=.(Biotype=biotype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(Biotype,Spliced)]
dt_x$RNA<-'Xist'
colnames(dt_x)[3]<-'Count'
dt_x$perct<-dt_x$Count*100/1000

dt_a<-as.data.table(top1000_airn)[,.N,by=.(Biotype=biotype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(Biotype,Spliced)]
dt_a$RNA<-'Airn'
colnames(dt_a)[3]<-'Count'
dt_a$perct<-dt_a$Count*100/1000

dt_k<-as.data.table(top1000_ot1)[,.N,by=.(Biotype=biotype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(Biotype,Spliced)]
dt_k$RNA<-'Kcnq1ot1'
colnames(dt_k)[3]<-'Count'
dt_k$perct<-dt_k$Count*100/1000

dt_c<-as.data.table(comb)[,.N,by=.(Biotype=biotype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(Biotype,Spliced)]
dt_c$RNA<-'Chrom'
colnames(dt_c)[3]<-'Count'
dt_c$perct<-dt_c$Count*100/19295

dt_comb<-rbind(dt_c,dt_x,dt_a,dt_k)

dt_comb$Biotype[which(dt_comb$Biotype=='protein_coding')]<-'protein\ncoding'
#Now want to graph percentages of the spliced and unspliced lncRNAs

dt_comb$RNA <- factor(dt_comb$RNA, levels = c("Chrom","Airn", "Kcnq1ot1", "Xist"))
dt_comb$Biotype <- factor(dt_comb$Biotype, levels = c("lncRNA","protein\ncoding", "others"))


library(LaCroixColoR)
pdf("XAK_top1000_RIPr_Biotype_splice.pdf", width = 10, height = 10)
ggplot(dt_comb) +
  geom_bar(aes(x = RNA, y = perct, fill = Spliced),
           position = "stack",
           stat = "identity") +
  facet_grid(~ Biotype, switch = "x") +
  labs(y = "Percent of Transcripts") +
  coord_cartesian(ylim=c(0,100))+
  theme(strip.placement = "outside",
        panel.spacing = unit(-.01,"cm"),
        plot.background = element_rect(fill = "transparent", color = NA),  # transparent background
        panel.background = element_rect(fill = "transparent", color = NA), # transparent panel background
        text = element_text(size = 36),  # Arial font with size 18
        strip.background = element_rect(fill = "transparent", color = NA), # transparent facet strips
        strip.text = element_text(size = 36),
        legend.title=element_blank(),
        legend.position = 'top',
        axis.title.x = element_blank(),
        axis.text.x=element_text(angle = 45, vjust=1,hjust=1)
        ) +
  scale_fill_manual(
    values = c("IntronExcluding" = "#F6A1A5", "IntronIncluding" = "#1BB6AF")  # Custom colors
  )


dev.off()



###############################
# Fisher exact test
# whether count of transcripts from lncRNA vs protein coding differ in the top 1000 for AKX vs in the total chromatin-associated population

#            pc vs lncRNA
# XAKtop1000
# chrom

typecount<-data.frame(RNA=c('Xist','Airn','Kcnq1ot1'),
                      pc_xak=as.numeric(joined_data[which(joined_data$Biotype=='protein_coding'),3:5]),
                      lnc_xak=as.numeric(joined_data[which(joined_data$Biotype=='lncRNA'),3:5]),
                      pc_chrom=joined_data[which(joined_data$Biotype=='protein_coding'),2],
                      lnc_chrom=joined_data[which(joined_data$Biotype=='lncRNA'),2])



typeFisher <- typecount %>% 
  crossing(comparison=c('greater','less')) %>%
  rowwise() %>%
  mutate(pval=fisher.test(
    matrix(c(pc_xak,lnc_xak,pc_chrom,lnc_chrom),
           nrow=2,byrow=T),
    alternative=comparison
  )$p.value) %>%
  ungroup()

typeFisher<-typeFisher%>%
  mutate(padj=p.adjust(pval,method='BH'))


write.csv(typeFisher,'XAK_RIPr_biotype_Fisher_stats.csv',row.names = F)


# within each biotype, whether the proportions of intron-excluding vs including differ in the top 1000 for AKX vs the total chromatin-associated population
dt_comb<-rbind(dt_c,dt_x,dt_a,dt_k)
dt_comb$RNA <- factor(dt_comb$RNA, levels = c("Chrom","Xist", "Airn", "Kcnq1ot1"))
dt_comb$Biotype <- factor(dt_comb$Biotype, levels = c("lncRNA","protein_coding", "others"))
dt_comb<-dt_comb[order(dt_comb$Biotype),]

splicecount<-data.frame(Biotype=rep(c('lncRNA','protein_coding','others'),each=3),
                        RNA=rep(c('Xist','Airn','Kcnq1ot1'),times=3),
                        IntInc_xak=as.numeric(dt_comb$Count[which(dt_comb$Spliced=='IntronIncluding' & dt_comb$RNA!='Chrom')]),
                        IntExc_xak=as.numeric(dt_comb$Count[which(dt_comb$Spliced=='IntronExcluding' & dt_comb$RNA!='Chrom')]),
                        IntInc_chrom=rep(dt_comb$Count[which(dt_comb$Spliced=='IntronIncluding' & dt_comb$RNA=='Chrom')],each=3),
                        IntExc_chrom=rep(dt_comb$Count[which(dt_comb$Spliced=='IntronExcluding' & dt_comb$RNA=='Chrom')],each=3)
                        )


spliceFisher <- splicecount %>% 
  crossing(comparison=c('greater','less')) %>%
  rowwise() %>%
  mutate(pval=fisher.test(
    matrix(c(IntInc_xak,IntExc_xak,IntInc_chrom,IntExc_chrom),
           nrow=2,byrow=T),
    alternative=comparison
  )$p.value) %>%
  ungroup()

spliceFisher<-spliceFisher%>%
  mutate(padj=p.adjust(pval,method='BH'))


write.csv(spliceFisher,'XAK_RIPr_splice_Fisher_stats.csv',row.names = F)



