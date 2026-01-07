# for CLIP CLAP
# get signal to noise bar plot for X and K 
# for HK and HU calculate rpmlessCK and CK and plot barplot
# sum all rpms under true peaks that falls in X and K

clap_hg<-read.csv('hg38CLAP_all_rpkm.csv',header=T)
clap_mm<-read.csv('mm10CLAP_10xk_rpkm.csv',header=T)

clip_hg<-read.csv('hg38CLAP10xk_CLIP_rpkm.csv',header=T)
clip_mm<-read.csv('mm10CLAP10xk_CLIP_rpkm.csv',header=T)


# human
# XIST chrX:73820656-73852723
# KOT1: chr11:2608328-2699994

get_hx<-function(df) {
  xidx<-which(df$Chr=='chrX' & df$Strand=='-' & df$Start>=73820656 & df$End<73852723)
  return(df[xidx,])
}

get_hot1<-function(df){
  kidx<-which(df$Chr=='chr11' & df$Strand=='-' & df$Start>=2608328 & df$End<2699994)
  return(df[kidx,])
}

# mouse
# Xist: chrX:103460366-103483254
# Kot1: chr7:143203458-143296549

get_mx<-function(df){
  xidx<-which(df$Chr=='chrX' & df$Strand=='-' & df$Start>=103460366 & df$End<103483254)
  return(df[xidx,])
}

get_mot1<-function(df){
  kidx<-which(df$Chr=='chr7' & df$Strand=='-' & df$Start>=143203458 & df$End<143296549)
  return(df[kidx,])
}


# get all stats for Xist
xstats<-matrix(nrow=12,ncol=4)
xstats<-as.data.frame(xstats)
colnames(xstats)<-c('condition','replicate','groups','rpmsum')
xstats$condition<-rep(c('CLIP_mm10_HNRNPU','CLIP_hg38_HNRNPU','CLAP_mm10_HNRNPU','CLAP_hg38_HNRNPU'),each=3)
xstats$replicate<-rep(c('r1','r1','r2','r1','r2','r1'),times=2)
xstats$groups<-rep(c('RBPlessCK','CK','CK','RBPlessCK','RBPlessCK','CK'),times=2)


clip_mmx<-get_mx(clip_mm)
clip_hgx<-get_hx(clip_hg)
clap_mmx<-get_mx(clap_mm)
clap_hgx<-get_hx(clap_hg)

  
xstats$rpmsum[1]<-sum(clip_mmx$mm10CLAP10k_CLIP_rpm_lessCK)
xstats$rpmsum[2]<-sum(clip_mmx$mm10CLAP10k_CLIP_ck1_rpm)
xstats$rpmsum[3]<-sum(clip_mmx$mm10CLAP10k_CLIP_ck2_rpm)

xstats$rpmsum[4]<-sum(clip_hgx$hg38CLAP10k_CLIP_Rep1_rpm_lessCK)
xstats$rpmsum[5]<-sum(clip_hgx$hg38CLAP10k_CLIP_Rep2_rpm_lessCK)
xstats$rpmsum[6]<-sum(clip_hgx$hg38CLAP10k_CLIP_ck_rpm)

xstats$rpmsum[7]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_rpm_lessCK)
xstats$rpmsum[8]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_ck1_rpm)
xstats$rpmsum[9]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_ck2_rpm)

xstats$rpmsum[10]<-sum(clap_hgx$SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK)
xstats$rpmsum[11]<-sum(clap_hgx$SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK)
xstats$rpmsum[12]<-sum(clap_hgx$SAFA_PlusTag_CLAP_hg38_ck_rpm)

write.csv(xstats,'Xist_CLIPCLAP_signaltonoise_points.csv',row.names = F)


# get all stats for Kot1
xstats<-matrix(nrow=12,ncol=4)
xstats<-as.data.frame(xstats)
colnames(xstats)<-c('condition','replicate','groups','rpmsum')
xstats$condition<-rep(c('CLIP_mm10_HNRNPU','CLIP_hg38_HNRNPU','CLAP_mm10_HNRNPU','CLAP_hg38_HNRNPU'),each=3)
xstats$replicate<-rep(c('r1','r1','r2','r1','r2','r1'),times=2)
xstats$groups<-rep(c('RBPlessCK','CK','CK','RBPlessCK','RBPlessCK','CK'),times=2)


clip_mmx<-get_mot1(clip_mm)
clip_hgx<-get_hot1(clip_hg)
clap_mmx<-get_mot1(clap_mm)
clap_hgx<-get_hot1(clap_hg)


xstats$rpmsum[1]<-sum(clip_mmx$mm10CLAP10k_CLIP_rpm_lessCK)
xstats$rpmsum[2]<-sum(clip_mmx$mm10CLAP10k_CLIP_ck1_rpm)
xstats$rpmsum[3]<-sum(clip_mmx$mm10CLAP10k_CLIP_ck2_rpm)

xstats$rpmsum[4]<-sum(clip_hgx$hg38CLAP10k_CLIP_Rep1_rpm_lessCK)
xstats$rpmsum[5]<-sum(clip_hgx$hg38CLAP10k_CLIP_Rep2_rpm_lessCK)
xstats$rpmsum[6]<-sum(clip_hgx$hg38CLAP10k_CLIP_ck_rpm)

xstats$rpmsum[7]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_rpm_lessCK)
xstats$rpmsum[8]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_ck1_rpm)
xstats$rpmsum[9]<-sum(clap_mmx$SAFA_MinusTag_CLAP_mm10_ck2_rpm)

xstats$rpmsum[10]<-sum(clap_hgx$SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK)
xstats$rpmsum[11]<-sum(clap_hgx$SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK)
xstats$rpmsum[12]<-sum(clap_hgx$SAFA_PlusTag_CLAP_hg38_ck_rpm)

write.csv(xstats,'Kot1_CLIPCLAP_signaltonoise_points.csv',row.names = F)


######### manually merge CLIPCLAP and RIP data together

library(ggplot2)
library(dplyr)

pratio<-read.csv('Xist_CLIPCLAPRIP_signaltonoise_means.csv',header=T)
ppoint<-read.csv('Xist_CLIPCLAPRIP_signaltonoise_points.csv',header=T)


pratio  <- pratio  %>% mutate(
  condition = factor(condition, levels=c('RIP_mm10_HNRNPK','RIP_mm10_HNRNPU',
                                         'CLIP_mm10_HNRNPU','CLIP_hg38_HNRNPU',
                                         'CLAP_mm10_HNRNPU','CLAP_hg38_HNRNPU')),
  groups = factor(groups,levels=c('RBPlessCK','CK'))
)
ppoint  <- ppoint  %>% mutate(
  condition = factor(condition, levels=levels(pratio$condition)),
  groups  = factor(groups,  levels = levels(pratio$groups))
)

## use a common dodging position so bars and points align
pos_bar <- position_dodge(width = 0.8)
pos_pts <- position_jitterdodge(dodge.width = 0.8, jitter.width = 0.35, jitter.height = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = pratio,
    aes(
      x = condition,
      y = rpmsum,
      fill  = groups,                           
      group = interaction(condition, groups)       # defines dodging groups
    ),
    position = pos_bar,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = ppoint,
    aes(
      x = condition,
      y = rpmsum,                          
      group = interaction(condition, groups),
      colour = groups# same grouping for same dodge
    ),
    position = pos_pts,                      # EXACTLY the same position object
    size     = 3.0,
    alpha    = 0.6
  ) +
  # some typical cosmetics
  labs(x = 'Samples', y = "RPM Sum") +
  scale_fill_manual(values = c('#4f4a8c','#9d9d9d'))+
  scale_color_manual(values = c('black','black'))+
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=5, r=5, b=5, l=60, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='right',
        legend.text=element_text(size=26),
        axis.text.x=element_text(size=26,color='black',angle=45,hjust=1),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))


ggsave('Xist_CLIPCLAPRIP_signaltonoise_barplot.pdf',plot=p,device='pdf',width=12,height=5.5,unit='in') 

#########################################
pratio<-read.csv('Kot1_CLIPCLAPRIP_signaltonoise_means.csv',header=T)
ppoint<-read.csv('Kot1_CLIPCLAPRIP_signaltonoise_points.csv',header=T)


pratio  <- pratio  %>% mutate(
  condition = factor(condition, levels=c('RIP_mm10_HNRNPK','RIP_mm10_HNRNPU',
                                         'CLIP_mm10_HNRNPU','CLIP_hg38_HNRNPU',
                                         'CLAP_mm10_HNRNPU','CLAP_hg38_HNRNPU')),
  groups = factor(groups,levels=c('RBPlessCK','CK'))
)
ppoint  <- ppoint  %>% mutate(
  condition = factor(condition, levels=levels(pratio$condition)),
  groups  = factor(groups,  levels = levels(pratio$groups))
)

## use a common dodging position so bars and points align
pos_bar <- position_dodge(width = 0.8)
pos_pts <- position_jitterdodge(dodge.width = 0.8, jitter.width = 0.35, jitter.height = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = pratio,
    aes(
      x = condition,
      y = rpmsum+1,
      fill  = groups,                           
      group = interaction(condition, groups)       # defines dodging groups
    ),
    position = pos_bar,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = ppoint,
    aes(
      x = condition,
      y = rpmsum+1,                          
      group = interaction(condition, groups),
      color = groups
    ),
    position = pos_pts,                      # EXACTLY the same position object
    size     = 3.0,
    alpha    = 0.6
  ) +
  # some typical cosmetics
  labs(x = 'Samples', y = "RPM Sum") +
  scale_fill_manual(values = c('#4f4a8c','#9d9d9d'))+
  scale_color_manual(values = c('black','black'))+
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=5, r=5, b=5, l=60, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='right',
        legend.text=element_text(size=26),
        axis.text.x=element_text(size=26,color='black',angle=45,hjust=1),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))


ggsave('Kot1_CLIPCLAPRIP_signaltonoise_barplot.pdf',plot=p,device='pdf',width=12,height=5.5,unit='in') 

