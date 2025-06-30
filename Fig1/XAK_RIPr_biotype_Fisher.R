#genetype analysis first on the top 1000 transcripts similar to X, K, A by RBP binding intensity

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

comb<-read.csv('allgenes_XAKRIPr_XAKnetworkr_common_sil_modu_genetype_06302025.csv',header=T)

comb<-comb[,c(1:4,14)]

# edit genetype
# comb$genetype[which(comb$genetype!='protein_coding' & comb$genetype!='lncRNA')]<-'others'

# get top1000 xak, exclude XAK from the list
top1000_xist <- comb[order(comb$Xist_RIP_r,decreasing = TRUE)[2:1001],]
top1000_airn <- comb[order(comb$Airn_RIP_r,decreasing = TRUE)[2:1001],]
top1000_ot1 <- comb[order(comb$Kcnq1ot1_RIP_r,decreasing = TRUE)[2:1001],]
  

genetype_counts_xist <- as.data.frame(table(top1000_xist$genetype))
genetype_counts_airn <- as.data.frame(table(top1000_airn$genetype))
genetype_counts_ot1 <- as.data.frame(table(top1000_ot1$genetype))

genetype_counts_chrom <- as.data.frame(table(comb$genetype))

# Rename the columns
colnames(genetype_counts_xist) <- c("genetype", "Count_xist")
colnames(genetype_counts_ot1) <- c("genetype", "Count_ot1")
colnames(genetype_counts_airn) <- c("genetype", "Count_airn")
colnames(genetype_counts_chrom) <- c("genetype", "Count_chrom")

joined_data <- full_join(genetype_counts_chrom, genetype_counts_xist, by = "genetype") %>%
  full_join(genetype_counts_airn, by = "genetype")%>%
  full_join(genetype_counts_ot1, by = "genetype")


#use the code below to generate spliced and unspliced info for each genetype
dt_x<-as.data.table(top1000_xist)[,.N,by=.(genetype=genetype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(genetype,Spliced)]
dt_x$RNA<-'Xist'
colnames(dt_x)[3]<-'Count'
dt_x$perct<-dt_x$Count*100/1000

dt_a<-as.data.table(top1000_airn)[,.N,by=.(genetype=genetype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(genetype,Spliced)]
dt_a$RNA<-'Airn'
colnames(dt_a)[3]<-'Count'
dt_a$perct<-dt_a$Count*100/1000

dt_k<-as.data.table(top1000_ot1)[,.N,by=.(genetype=genetype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(genetype,Spliced)]
dt_k$RNA<-'Kcnq1ot1'
colnames(dt_k)[3]<-'Count'
dt_k$perct<-dt_k$Count*100/1000

dt_c<-as.data.table(comb)[,.N,by=.(genetype=genetype,Spliced=ifelse(grepl("unspliced",gene),'IntronIncluding','IntronExcluding'))][order(genetype,Spliced)]
dt_c$RNA<-'Chrom'
colnames(dt_c)[3]<-'Count'
dt_c$perct<-dt_c$Count*100/19295

dt_comb<-rbind(dt_c,dt_x,dt_a,dt_k)

dt_comb$genetype[which(dt_comb$genetype=='protein_coding')]<-'protein\ncoding'
#Now want to graph percentages of the spliced and unspliced lncRNAs

dt_comb$RNA <- factor(dt_comb$RNA, levels = c("Chrom","Airn", "Kcnq1ot1", "Xist"))
dt_comb$genetype <- factor(dt_comb$genetype, levels = c("lncRNA","protein\ncoding", "others"))


library(LaCroixColoR)
pdf("XAK_top1000_RIPr_genetype_splice.pdf", width = 10, height = 10)
ggplot(dt_comb) +
  geom_bar(aes(x = RNA, y = perct, fill = Spliced),
           position = "stack",
           stat = "identity") +
  facet_grid(~ genetype, switch = "x") +
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
                      pc_xak=as.numeric(joined_data[which(joined_data$genetype=='protein_coding'),3:5]),
                      lnc_xak=as.numeric(joined_data[which(joined_data$genetype=='lncRNA'),3:5]),
                      pc_chrom=joined_data[which(joined_data$genetype=='protein_coding'),2],
                      lnc_chrom=joined_data[which(joined_data$genetype=='lncRNA'),2])



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


write.csv(typeFisher,'XAK_RIPr_genetype_Fisher_stats.csv',row.names = F)


# within each genetype, whether the proportions of intron-excluding vs including differ in the top 1000 for AKX vs the total chromatin-associated population
dt_comb<-rbind(dt_c,dt_x,dt_a,dt_k)
dt_comb$RNA <- factor(dt_comb$RNA, levels = c("Chrom","Xist", "Airn", "Kcnq1ot1"))
dt_comb$genetype <- factor(dt_comb$genetype, levels = c("lncRNA","protein_coding", "others"))
dt_comb<-dt_comb[order(dt_comb$genetype),]

splicecount<-data.frame(genetype=rep(c('lncRNA','protein_coding','others'),each=3),
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



