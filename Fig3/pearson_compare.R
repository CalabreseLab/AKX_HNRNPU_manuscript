#compare inter vs intra community pearson distrubtion 

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)


############### X
dfx<-read.csv('Xist_allpeaks200_pearson.csv',header=T,row.names = 1)

combx <- dfx %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval") %>%
  # keep only those where the original row-index < col-index
  filter(
    match(peak1, rownames(dfx)) < match(peak2, names(dfx))
  )


combx$com1<-gsub('_p.*$','',combx$peak1)
combx$com2<-gsub('_p.*$','',combx$peak2)


table(combx$com1,combx$com2)

#manually remove Xist c7 as it only has one peak
combx<-combx[which(combx$com1!='Xist_c7' & combx$com2!='Xist_c7'),]

########## A
dfa<-read.csv('Airn_allpeaks200_pearson.csv',header=T,row.names = 1)

comba <- dfa %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval") %>%
  # keep only those where the original row-index < col-index
  filter(
    match(peak1, rownames(dfa)) < match(peak2, names(dfa))
  )


comba$com1<-gsub('_p.*$','',comba$peak1)
comba$com2<-gsub('_p.*$','',comba$peak2)

table(comba$com1,comba$com2)

# manually remove Airn c5 as it only has one minor peak
# comba<-comba[which(comba$com1!='Airn_c5' & comba$com2!='Airn_c5'),]

############ K

dfk<-read.csv('Kot1_allpeaks200_pearson.csv',header=T,row.names = 1)

combk <- dfk %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval") %>%
  # keep only those where the original row-index < col-index
  filter(
    match(peak1, rownames(dfk)) < match(peak2, names(dfk))
  )


combk$com1<-gsub('_p.*$','',combk$peak1)
combk$com2<-gsub('_p.*$','',combk$peak2)


table(combk$com1,combk$com2)


##################

comb<-rbind(combx,comba,combk)

comlist<-unique(comb$com1)

comstats<-data.frame(community=comlist,
                     intrasize=0,
                     intersize=0,
                     pval=NA)

complot<-data.frame(community=as.character(),
                    group=as.character(),
                    rval=as.numeric())

for (n in 1:length(comlist)) {
  
  cp<-comlist[n]
  
  intra<-comb$rval[which(comb$com1==cp & comb$com2==cp)]
  
  inter<-comb$rval[which((comb$com1==cp & comb$com2!=cp) | (comb$com1!=cp & comb$com2==cp))]
  
  comstats$intrasize[n]<-length(intra)
  comstats$intersize[n]<-length(inter)
  
  comstats$pval[n]<-wilcox.test(intra,inter,alternative = 'greater')$p.value
  
  
  tempintra<-data.frame(community=cp,group='intra',rval=intra)
  tempinter<-data.frame(community=cp,group='inter',rval=inter)
  
  complot<-rbind(complot,tempintra,tempinter)
  
}

write.csv(comstats,'allpeak_intrainter_wilcoxp.csv',row.names = F)

library(ggplot2)
library(gtools)

plotx<-mixedsort(comlist)
plotx<-plotx[c(17:21,1:16)]
complot$community<-factor(complot$community,levels=plotx)
complot$group<-factor(complot$group,levels=c('intra','inter'))

p<-ggplot(data=complot,aes(x=community,y=rval,fill=group,color=group)) +
  labs(x = "Community",y = "Pearson's r") +
  coord_cartesian(ylim=c(-0.4,0.6))+
  scale_y_continuous(breaks=seq(from=-0.4, to=0.6, by=0.2)) +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='top',
        legend.title=element_blank(),
        legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x=element_text(size=30, angle=90,hjust=0.5,vjust=0.5),
        axis.text.y=element_text(size=30)) +
  geom_boxplot(color='black',outlier.shape = NA,width = 0.6,position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')

ggsave("allpeak_intrainter_boxplot.pdf", plot = p, width = 18, height = 8)

###########################################
# compare to null

# get all the intra r vals from comb
combintra<-comb[which(comb$com1==comb$com2),]


##### X
nullx<-read.csv('Xist_peakvsnull_pearson.csv',header=T,row.names = 1)

nullx <- nullx %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval")
  # keep only those where the original row-index < col-index
  

nullx$com<-gsub('_p.*$','',nullx$peak1)


table(nullx$com)

#manually remove Xist c7 as it only has one peak
nullx<-nullx[which(nullx$com!='Xist_c7'),]

####### A
nulla<-read.csv('Airn_peakvsnull_pearson.csv',header=T,row.names = 1)

nulla <- nulla %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval")
# keep only those where the original row-index < col-index


nulla$com<-gsub('_p.*$','',nulla$peak1)


table(nulla$com)

# manually remove Airn c5 as it only has one minor peak
#nulla<-nulla[which(nulla$com!='Airn_c5'),]

####### K
nullk<-read.csv('Kot1_peakvsnull_pearson.csv',header=T,row.names = 1)

nullk <- nullk %>%
  # move rownames into a column
  rownames_to_column(var = "peak1") %>%
  # pivot to long form
  pivot_longer(-peak1, names_to = "peak2", values_to = "rval")
# keep only those where the original row-index < col-index


nullk$com<-gsub('_p.*$','',nullk$peak1)


table(nullk$com)

##########################################

nullall<-rbind(nullx,nulla,nullk)

comlist<-unique(combintra$com1)

comstats<-data.frame(community=comlist,
                     intrasize=0,
                     nullsize=0,
                     pval=NA)

complot<-data.frame(community=as.character(),
                    group=as.character(),
                    rval=as.numeric())

for (n in 1:length(comlist)) {
  
  cp<-comlist[n]
  
  intra<-combintra$rval[which(combintra$com1==cp & combintra$com2==cp)]
  
  nullr<-nullall$rval[which(nullall$com==cp)]
  
  comstats$intrasize[n]<-length(intra)
  comstats$nullsize[n]<-length(nullr)
  
  comstats$pval[n]<-wilcox.test(intra,nullr,alternative = 'greater')$p.value
  
  
  tempintra<-data.frame(community=cp,group='intra',rval=intra)
  tempnull<-data.frame(community=cp,group='null',rval=nullr)
  
  complot<-rbind(complot,tempintra,tempnull)
  
}

write.csv(comstats,'allpeak_intranull_wilcoxp.csv',row.names = F)

library(ggplot2)
library(gtools)

plotx<-mixedsort(comlist)
plotx<-plotx[c(17:21,1:16)]
complot$community<-factor(complot$community,levels=plotx)
complot$group<-factor(complot$group,levels=c('intra','null'))

p<-ggplot(data=complot,aes(x=community,y=rval,fill=group,color=group)) +
  labs(x = "Community",y = "Pearson's r") +
  coord_cartesian(ylim=c(-0.4,0.6))+
  scale_y_continuous(breaks=seq(from=-0.4, to=0.6, by=0.2)) +
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='top',
        legend.title=element_blank(),
        legend.text = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x=element_text(size=30, angle=90,hjust=0.5,vjust=0.5),
        axis.text.y=element_text(size=30)) +
  geom_boxplot(color='black',outlier.shape = NA,width = 0.6,position = position_dodge(width = 0.8)) +
  scale_fill_brewer(palette = 'Set1')+
  scale_color_brewer(palette = 'Set1')

ggsave("allpeak_intranull_boxplot.pdf", plot = p, width = 18, height = 8)

