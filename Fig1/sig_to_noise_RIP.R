
####################################################
# for hm_ samples
# get signal to noise box plot for top 10K 5K 2.5K and 1K
# for HK and HU calculate rpmlessCK and CK and plot barplot

mm10_hk<-read.csv('mm10_HNRNPK_all_rpkm.csv',header=T)
mm10_hu<-read.csv('mm10_HNRNPU_all_rpkm.csv',header=T)

######### HK

hk_mm_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(hk_mm_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
hk_mm_comb$exp<-'HNRNPK'
hk_mm_comb$group<-c('10K','5K','2.5K','1K')
hk_mm_comb$group<-factor(hk_mm_comb$group,levels=c('10K','5K','2.5K','1K'))
hk_mm_comb$sp<-'mm10'

temp<-mm10_hk[which(mm10_hk$ranking<=10000),c(12,13,21,22)]
tempsum<-colSums(temp)
hk_mm_comb[1,4]<-tempsum['hm_HNRNPK_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r1_rpm']
hk_mm_comb[1,5]<-tempsum['hm_HNRNPK_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r2_rpm']

temp<-mm10_hk[which(mm10_hk$ranking<=5000),c(12,13,21,22)]
tempsum<-colSums(temp)
hk_mm_comb[2,4]<-tempsum['hm_HNRNPK_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r1_rpm']
hk_mm_comb[2,5]<-tempsum['hm_HNRNPK_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r2_rpm']

temp<-mm10_hk[which(mm10_hk$ranking<=2500),c(12,13,21,22)]
tempsum<-colSums(temp)
hk_mm_comb[3,4]<-tempsum['hm_HNRNPK_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r1_rpm']
hk_mm_comb[3,5]<-tempsum['hm_HNRNPK_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r2_rpm']

temp<-mm10_hk[which(mm10_hk$ranking<=1000),c(12,13,21,22)]
tempsum<-colSums(temp)
hk_mm_comb[4,4]<-tempsum['hm_HNRNPK_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r1_rpm']
hk_mm_comb[4,5]<-tempsum['hm_HNRNPK_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPK_flag_r2_rpm']


######### HU

hu_mm_comb<-as.data.frame(matrix(nrow=4,ncol=5))
colnames(hu_mm_comb)<-c('exp','group','sp','stn_perct_r1','stn_perct_r2')
hu_mm_comb$exp<-'HNRNPU'
hu_mm_comb$group<-c('10K','5K','2.5K','1K')
hu_mm_comb$group<-factor(hu_mm_comb$group,levels=c('10K','5K','2.5K','1K'))
hu_mm_comb$sp<-'mm10'

temp<-mm10_hu[which(mm10_hu$ranking<=10000),c(12,13,21,22)]
tempsum<-colSums(temp)
hu_mm_comb[1,4]<-tempsum['hm_HNRNPU_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r1_rpm']
hu_mm_comb[1,5]<-tempsum['hm_HNRNPU_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r2_rpm']

temp<-mm10_hu[which(mm10_hu$ranking<=5000),c(12,13,21,22)]
tempsum<-colSums(temp)
hu_mm_comb[2,4]<-tempsum['hm_HNRNPU_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r1_rpm']
hu_mm_comb[2,5]<-tempsum['hm_HNRNPU_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r2_rpm']

temp<-mm10_hu[which(mm10_hu$ranking<=2500),c(12,13,21,22)]
tempsum<-colSums(temp)
hu_mm_comb[3,4]<-tempsum['hm_HNRNPU_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r1_rpm']
hu_mm_comb[3,5]<-tempsum['hm_HNRNPU_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r2_rpm']

temp<-mm10_hu[which(mm10_hu$ranking<=1000),c(12,13,21,22)]
tempsum<-colSums(temp)
hu_mm_comb[4,4]<-tempsum['hm_HNRNPU_flag_r1_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r1_rpm']
hu_mm_comb[4,5]<-tempsum['hm_HNRNPU_flag_r2_rpm_lessCK']*100/tempsum['hm_HNRNPU_flag_r2_rpm']

ripcomb<-rbind(hk_mm_comb,hu_mm_comb)

write.csv(ripcomb,'RIP_signaltonoise_raw.csv',row.names = F)

# manualy edits to make the points and ratios file
################################

library(ggplot2)
library(dplyr)

pratio<-read.csv('RIP_signaltonoise_ratios.csv',header=T)
ppoint<-read.csv('RIP_signaltonoise_points.csv',header=T)

## optional: make sure group/method are factors in a sensible order
pratio  <- pratio  %>% mutate(
  exp  = factor(exp,levels=c('HNRNPK','HNRNPU')),
  group = factor(group,levels=c('10K','5K','2.5K','1K'))
)
ppoint  <- ppoint  %>% mutate(
  exp  = factor(exp,  levels = levels(pratio$exp)),
  group = factor(group, levels = levels(pratio$group))
)

## use a common dodging position so bars and points align
pos_dodge <- position_dodge2(width = 0.8,preserve = 'single',padding = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = pratio,
    aes(
      x = exp,
      y = stn_perct_r12mean,
      fill  = group,                           
      group = interaction(group, exp)       # defines dodging groups
    ),
    position = pos_dodge,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = ppoint,
    aes(
      x = exp,
      y = stn_perct, 
      group = interaction(group,exp)       
    ),
    color = 'black',
    position = pos_dodge,
    size     = 3.0,
    alpha    = 0.6
  ) +
  labs(x = 'RIP', y = "SignaltoNoise%") +
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


ggsave('RIP_signaltonoise_barplot.pdf',plot=p,device='pdf',width=4.5,height=4.5,unit='in') 
