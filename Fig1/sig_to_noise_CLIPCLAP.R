
#####################################################
# for CLIP CLAP
# get signal to noise boxplot for top 10k 5k 2.5k and 1k
# for human and mouse calculate rpmlessCK and CK and plot barplot
# sum all rpms under true peaks 

library(tidyverse)
library(ggplot2)


#############################################
clap_hg<-read.csv('hg38CLAP_all_rpkm.csv',header=T)
clap_mm<-read.csv('mm10CLAP_10xk_rpkm.csv',header=T)


# human
clap_hg_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(clap_hg_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
clap_hg_comb$exp<-'CLAP'
clap_hg_comb$group<-c('10K','5K','2.5K','1K')
clap_hg_comb$group<-factor(clap_hg_comb$group,levels=c('10K','5K','2.5K','1K'))
clap_hg_comb$sp<-'hg38'


temp<-clap_hg[which(clap_hg$ranking<=10000),c(12,13,23,24)]
tempsum<-colSums(temp)
clap_hg_comb[1,4]<-tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm']
clap_hg_comb[1,5]<-tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm']

temp<-clap_hg[which(clap_hg$ranking<=5000),c(12,13,23,24)]
tempsum<-colSums(temp)
clap_hg_comb[2,4]<-tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm']
clap_hg_comb[2,5]<-tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm']

temp<-clap_hg[which(clap_hg$ranking<=2500),c(12,13,23,24)]
tempsum<-colSums(temp)
clap_hg_comb[3,4]<-tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm']
clap_hg_comb[3,5]<-tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm']

temp<-clap_hg[which(clap_hg$ranking<=1000),c(12,13,23,24)]
tempsum<-colSums(temp)
clap_hg_comb[4,4]<-tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep1_hg38_rpm']
clap_hg_comb[4,5]<-tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK']*100/tempsum['SAFA_PlusTag_CLAP_Rep2_hg38_rpm']


##### mouse

clap_mm_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(clap_mm_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
clap_mm_comb$exp<-'CLAP'
clap_mm_comb$group<-c('10K','5K','2.5K','1K')
clap_mm_comb$group<-factor(clap_mm_comb$group,levels=c('10K','5K','2.5K','1K'))
clap_mm_comb$sp<-'mm10'


temp<-clap_mm[which(clap_mm$ranking<=10000),c(9,11)]
tempsum<-colSums(temp)
clap_mm_comb[1,4]<-tempsum['SAFA_MinusTag_CLAP_mm10_rpm_lessCK']*100/tempsum['SAFA_MinusTag_CLAP_mm10_rpm']

temp<-clap_mm[which(clap_mm$ranking<=5000),c(9,11)]
tempsum<-colSums(temp)
clap_mm_comb[2,4]<-tempsum['SAFA_MinusTag_CLAP_mm10_rpm_lessCK']*100/tempsum['SAFA_MinusTag_CLAP_mm10_rpm']

temp<-clap_mm[which(clap_mm$ranking<=2500),c(9,11)]
tempsum<-colSums(temp)
clap_mm_comb[3,4]<-tempsum['SAFA_MinusTag_CLAP_mm10_rpm_lessCK']*100/tempsum['SAFA_MinusTag_CLAP_mm10_rpm']

temp<-clap_mm[which(clap_mm$ranking<=1000),c(9,11)]
tempsum<-colSums(temp)
clap_mm_comb[4,4]<-tempsum['SAFA_MinusTag_CLAP_mm10_rpm_lessCK']*100/tempsum['SAFA_MinusTag_CLAP_mm10_rpm']

clap_comb<-rbind(clap_hg_comb,clap_mm_comb)


#################################################
clip_hg<-read.csv('hg38CLAP10xk_CLIP_rpkm.csv',header=T)
clip_mm<-read.csv('mm10CLAP10xk_CLIP_rpkm.csv',header=T)

# human
clip_hg_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(clip_hg_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
clip_hg_comb$exp<-'CLIP'
clip_hg_comb$group<-c('10K','5K','2.5K','1K')
clip_hg_comb$group<-factor(clip_hg_comb$group,levels=c('10K','5K','2.5K','1K'))
clip_hg_comb$sp<-'hg38'


temp<-clip_hg[which(clip_hg$Geneid<=10000),c(12,13,16,17)]
tempsum<-colSums(temp)
clip_hg_comb[1,4]<-tempsum['hg38CLAP10k_CLIP_Rep1_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep1_rpm']
clip_hg_comb[1,5]<-tempsum['hg38CLAP10k_CLIP_Rep2_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep2_rpm']

temp<-clip_hg[which(clip_hg$Geneid<=5000),c(12,13,16,17)]
tempsum<-colSums(temp)
clip_hg_comb[2,4]<-tempsum['hg38CLAP10k_CLIP_Rep1_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep1_rpm']
clip_hg_comb[2,5]<-tempsum['hg38CLAP10k_CLIP_Rep2_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep2_rpm']

temp<-clip_hg[which(clip_hg$Geneid<=2500),c(12,13,16,17)]
tempsum<-colSums(temp)
clip_hg_comb[3,4]<-tempsum['hg38CLAP10k_CLIP_Rep1_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep1_rpm']
clip_hg_comb[3,5]<-tempsum['hg38CLAP10k_CLIP_Rep2_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep2_rpm']

temp<-clip_hg[which(clip_hg$Geneid<=1000),c(12,13,16,17)]
tempsum<-colSums(temp)
clip_hg_comb[4,4]<-tempsum['hg38CLAP10k_CLIP_Rep1_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep1_rpm']
clip_hg_comb[4,5]<-tempsum['hg38CLAP10k_CLIP_Rep2_rpm_lessCK']*100/tempsum['hg38CLAP10k_CLIP_Rep2_rpm']


##### mouse

clip_mm_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(clip_mm_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
clip_mm_comb$exp<-'CLIP'
clip_mm_comb$group<-c('10K','5K','2.5K','1K')
clip_mm_comb$group<-factor(clip_mm_comb$group,levels=c('10K','5K','2.5K','1K'))
clip_mm_comb$sp<-'mm10'


temp<-clip_mm[which(clip_mm$Geneid<=10000),c(11,13)]
tempsum<-colSums(temp)
clip_mm_comb[1,4]<-tempsum['mm10CLAP10k_CLIP_rpm_lessCK']*100/tempsum['mm10CLAP10k_CLIP_rpm']

temp<-clip_mm[which(clip_mm$Geneid<=5000),c(11,13)]
tempsum<-colSums(temp)
clip_mm_comb[2,4]<-tempsum['mm10CLAP10k_CLIP_rpm_lessCK']*100/tempsum['mm10CLAP10k_CLIP_rpm']

temp<-clip_mm[which(clip_mm$Geneid<=2500),c(11,13)]
tempsum<-colSums(temp)
clip_mm_comb[3,4]<-tempsum['mm10CLAP10k_CLIP_rpm_lessCK']*100/tempsum['mm10CLAP10k_CLIP_rpm']

temp<-clip_mm[which(clip_mm$Geneid<=1000),c(11,13)]
tempsum<-colSums(temp)
clip_mm_comb[4,4]<-tempsum['mm10CLAP10k_CLIP_rpm_lessCK']*100/tempsum['mm10CLAP10k_CLIP_rpm']

clip_comb<-rbind(clip_hg_comb,clip_mm_comb)

ccomb<-rbind(clap_comb,clip_comb)

write.csv(ccomb,'CLIPCLAP_signaltonoise_raw.csv',row.names = F)

# manualy edits to make the points and ratios file
################################

library(ggplot2)
library(dplyr)

pratio<-read.csv('CLIPCLAP_signaltonoise_ratios.csv',header=T)
ppoint<-read.csv('CLIPCLAP_signaltonoise_points.csv',header=T)

## optional: make sure group/method are factors in a sensible order
pratio  <- pratio  %>% mutate(
  sp = factor(sp, levels=c('hg38','mm10')),
  exp  = factor(exp,levels=c('CLAP','CLIP')),
  group = factor(group,levels=c('10K','5K','2.5K','1K'))
)
ppoint  <- ppoint  %>% mutate(
  sp = factor(sp, levels=levels(pratio$sp)),
  exp  = factor(exp,  levels = levels(pratio$exp)),
  group = factor(group, levels = levels(pratio$group))
)

subpratio<-pratio[which(pratio$exp=='CLIP'),]
subppoint<-ppoint[which(ppoint$exp=='CLIP'),]

## use a common dodging position so bars and points align
pos_dodge <- position_dodge2(width = 0.8,preserve = 'single',padding = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = subpratio,
    aes(
      x = sp,
      y = stn_perct_r12mean,
      fill  = group,                           
      group = interaction(group, sp)       # defines dodging groups
    ),
    position = pos_dodge,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = subppoint,
    aes(
      x = sp,
      y = stn_perct,                          
      group = interaction(group,sp)       # same grouping for same dodge
    ),
    color = 'black',
    position = pos_dodge,                      # EXACTLY the same position object
    size     = 3.0,
    alpha    = 0.6
  ) +
  # some typical cosmetics
  labs(x = 'CLIP', y = "SignaltoNoise%") +
  coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = c('#393564','#4f4a8c','#655faa','#8581bc'))+
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=5, r=5, b=5, l=5, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='top',
        legend.text=element_text(size=26),
        axis.text.x=element_text(size=26,color='black'),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))


ggsave('CLIP_signaltonoise_barplot.pdf',plot=p,device='pdf',width=4.5,height=4.5,unit='in') 



subpratio<-pratio[which(pratio$exp=='CLAP'),]
subppoint<-ppoint[which(ppoint$exp=='CLAP'),]

## use a common dodging position so bars and points align
pos_dodge <- position_dodge2(width = 0.8,preserve = 'single',padding = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = subpratio,
    aes(
      x = sp,
      y = stn_perct_r12mean,
      fill  = group,                           
      group = interaction(group, sp)       # defines dodging groups
    ),
    position = pos_dodge,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = subppoint,
    aes(
      x = sp,
      y = stn_perct,                          
      group = interaction(group,sp)       # same grouping for same dodge
    ),
    color = 'black',
    position = pos_dodge,                      # EXACTLY the same position object
    size     = 3.0,
    alpha    = 0.6
  ) +
  # some typical cosmetics
  labs(x = 'CLAP', y = "SignaltoNoise%") +
  coord_cartesian(ylim = c(0,100))+
  scale_fill_manual(values = c('#393564','#4f4a8c','#655faa','#8581bc'))+
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=5, r=5, b=5, l=5, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='top',
        legend.text=element_text(size=26),
        axis.text.x=element_text(size=26,color='black'),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))


ggsave('CLAP_signaltonoise_barplot.pdf',plot=p,device='pdf',width=4.5,height=4.5,unit='in') 
