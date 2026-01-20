# organize kallisto data and calculate rpm-ck and ckrpm of each transcripts
# For paired-end data, kallisto reports fragment counts, i.e. read pairs, not individual reads

#################### human

#### CLAP
claphg<-read.table('SAFA_PlusTag_CLAP_hg38_peaks_kallisto/abundance.tsv',sep='\t',header=T)
claphgck<-read.table('SAFA_PlusTag_CLAP_hg38_CKpeaks_kallisto/abundance.tsv',sep='\t',header=T)

claphg$tpm<-NULL
claphgck<-claphgck[,c(1,4)]

colnames(claphg)[4]<-'CLAPhg_est_counts'
colnames(claphgck)[2]<-'CLAPhgCK_est_counts'

sum(claphg$target_id!=claphgck$target_id)
claphgcomb<-claphg
claphgcomb$CLAPhgCK_est_counts<-claphgck$CLAPhgCK_est_counts

# total fastq counts come from previous count file
claphgcomb$CLAPhg_rpm<-claphgcomb$CLAPhg_est_counts*1000000/44467024
claphgcomb$CLAPhgCK_rpm<-claphgcomb$CLAPhgCK_est_counts*1000000/5703926

claphgcomb$CLAPhg_rpm_lessControl<-claphgcomb$CLAPhg_rpm-claphgcomb$CLAPhgCK_rpm
claphgcomb$CLAPhg_rpm_lessControl[which(claphgcomb$CLAPhg_rpm_lessControl<0)]<-0

########### CLIP
cliphg<-read.table('SAFA_PlusTag_CLIP_hg38_peaks_kallisto/abundance.tsv',sep='\t',header=T)
cliphgck<-read.table('SAFA_PlusTag_CLIP_hg38_CKpeaks_kallisto/abundance.tsv',sep='\t',header=T)

cliphg$tpm<-NULL
cliphgck<-cliphgck[,c(1,4)]

colnames(cliphg)[4]<-'CLIPhg_est_counts'
colnames(cliphgck)[2]<-'CLIPhgCK_est_counts'

sum(cliphg$target_id!=cliphgck$target_id)
cliphgcomb<-cliphg
cliphgcomb$CLIPhgCK_est_counts<-cliphgck$CLIPhgCK_est_counts

cliphgcomb$CLIPhg_rpm<-cliphgcomb$CLIPhg_est_counts*1000000/38511502
cliphgcomb$CLIPhgCK_rpm<-cliphgcomb$CLIPhgCK_est_counts*1000000/4879320

cliphgcomb$CLIPhg_rpm_lessControl<-cliphgcomb$CLIPhg_rpm-cliphgcomb$CLIPhgCK_rpm
cliphgcomb$CLIPhg_rpm_lessControl[which(cliphgcomb$CLIPhg_rpm_lessControl<0)]<-0

###########################################
sum(claphgcomb$target_id!=cliphgcomb$target_id)
hgcomb<-cbind(claphgcomb,cliphgcomb[,c(4:8)])

xkhg<-hgcomb[which(hgcomb$target_id=='XIST_chrX_73820649_73852723_-_ENSG00000229807.14_ENST00000429829.6' | 
                     hgcomb$target_id=='KCNQ1OT1_chr11_2608328_2699994_-_ENSG00000269821.2_ENST00000597346.1.monoexonic.unspliced'),]

write.csv(xkhg,'XK_hg38.csv',row.names = F)


################################ mouse ################

##### CLAP
clapmm<-read.table('SAFA_MinusTag_CLAP_mm10_peaks_kallisto/abundance.tsv',sep='\t',header=T)
clapmmck<-read.table('SAFA_MinusTag_CLAP_mm10_CKpeaks_kallisto/abundance.tsv',sep='\t',header=T)

clapmm$tpm<-NULL
clapmmck<-clapmmck[,c(1,4)]

colnames(clapmm)[4]<-'CLAPmm_est_counts'
colnames(clapmmck)[2]<-'CLAPmmCK_est_counts'

sum(clapmm$target_id!=clapmmck$target_id)

clapmm$CLAPmmCK_est_counts<-clapmmck$CLAPmmCK_est_counts

clapmm$CLAPmm_rpm<-clapmm$CLAPmm_est_counts*1000000/5703926
clapmm$CLAPmmCK_rpm<-clapmm$CLAPmmCK_est_counts*1000000/44467024

clapmm$CLAPmm_rpm_lessControl<-clapmm$CLAPmm_rpm-clapmm$CLAPmmCK_rpm
clapmm$CLAPmm_rpm_lessControl[which(clapmm$CLAPmm_rpm_lessControl<0)]<-0

##### CLIP
clipmm<-read.table('SAFA_MinusTag_CLIP_mm10_peaks_kallisto/abundance.tsv',sep='\t',header=T)
clipmmck<-read.table('SAFA_MinusTag_CLIP_mm10_CKpeaks_kallisto/abundance.tsv',sep='\t',header=T)

clipmm$tpm<-NULL
clipmmck<-clipmmck[,c(1,4)]

colnames(clipmm)[4]<-'CLIPmm_est_counts'
colnames(clipmmck)[2]<-'CLIPmmCK_est_counts'

sum(clipmm$target_id!=clipmmck$target_id)

clipmm$CLIPmmCK_est_counts<-clipmmck$CLIPmmCK_est_counts

clipmm$CLIPmm_rpm<-clipmm$CLIPmm_est_counts*1000000/4879320
clipmm$CLIPmmCK_rpm<-clipmm$CLIPmmCK_est_counts*1000000/38511502

clipmm$CLIPmm_rpm_lessControl<-clipmm$CLIPmm_rpm-clipmm$CLIPmmCK_rpm
clipmm$CLIPmm_rpm_lessControl[which(clipmm$CLIPmm_rpm_lessControl<0)]<-0

###### RIP HK
hk<-read.table('hm_HNRNPK_kallisto/abundance.tsv',sep='\t',header=T)
hkck<-read.table('hm_HNRNPK_CK_kallisto/abundance.tsv',sep='\t',header=T)

hk$tpm<-NULL
hkck<-hkck[,c(1,4)]

colnames(hk)[4]<-'RIPHNRNPK_est_counts'
colnames(hkck)[2]<-'RIPHNRNPKCK_est_counts'

sum(hk$target_id!=hkck$target_id)

hk$RIPHNRNPKCK_est_counts<-hkck$RIPHNRNPKCK_est_counts

hk$RIPHNRNPK_rpm<-hk$RIPHNRNPK_est_counts*1000000/62978874
hk$RIPHNRNPKCK_rpm<-hk$RIPHNRNPKCK_est_counts*1000000/60319397

hk$RIPHNRNPK_rpm_lessControl<-hk$RIPHNRNPK_rpm-hk$RIPHNRNPKCK_rpm
hk$RIPHNRNPK_rpm_lessControl[which(hk$RIPHNRNPK_rpm_lessControl<0)]<-0

###### RIP HU
hu<-read.table('hm_HNRNPU_kallisto/abundance.tsv',sep='\t',header=T)
huck<-read.table('hm_HNRNPU_CK_kallisto/abundance.tsv',sep='\t',header=T)

hu$tpm<-NULL
huck<-huck[,c(1,4)]

colnames(hu)[4]<-'RIPHNRNPU_est_counts'
colnames(huck)[2]<-'RIPHNRNPUCK_est_counts'

sum(hu$target_id!=huck$target_id)

hu$RIPHNRNPUCK_est_counts<-huck$RIPHNRNPUCK_est_counts

hu$RIPHNRNPU_rpm<-hu$RIPHNRNPU_est_counts*1000000/69381631
hu$RIPHNRNPUCK_rpm<-hu$RIPHNRNPUCK_est_counts*1000000/60319397

hu$RIPHNRNPU_rpm_lessControl<-hu$RIPHNRNPU_rpm-hu$RIPHNRNPUCK_rpm
hu$RIPHNRNPU_rpm_lessControl[which(hu$RIPHNRNPU_rpm_lessControl<0)]<-0

############
sum(clapmm$target_id!=clipmm$target_id)
sum(clapmm$target_id!=hk$target_id)
sum(clapmm$target_id!=hu$target_id)
mmcomb<-cbind(clapmm,clipmm[,c(4:10)],hk[,c(4:10)],hu[,c(4:10)])


xkmm<-mmcomb[which(mmcomb$target_id=='Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3' | 
                     mmcomb$target_id=='Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced'),]

write.csv(xkmm,'XK_mm10.csv',row.names = F)


######################################################################
# manually organize the data to integrate XK hg38 and mm10 together
# manually merge CLIPCLAP and RIP data together

library(ggplot2)
library(dplyr)

pratio<-read.csv('Xist_CLIPCLAPRIP_signaltonoise_kallisto_means.csv',header=T)

pratio  <- pratio  %>% mutate(
  condition = factor(condition, levels=c('RIP_mm10_HNRNPK','RIP_mm10_HNRNPU',
                                         'CLIP_hg38_HNRNPU','CLIP_mm10_HNRNPU',
                                         'CLAP_hg38_HNRNPU','CLAP_mm10_HNRNPU')),
  groups = factor(groups,levels=c('RBPlessCK','CK'))
)


pos_bar <- position_dodge(width = 0.8)

p<-ggplot() +
  geom_col(
    data = pratio,
    aes(
      x = condition,
      y = rpmsum,
      fill  = groups,                           
      group = interaction(condition, groups)       
    ),
    position = pos_bar,
    width    = 0.8,
    linewidth = 1.0
  ) +
  labs(x = 'Samples', y = "RPM Sum") +
  scale_fill_manual(values = c('#4f4a8c','#9d9d9d'))+
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


ggsave('Xist_CLIPCLAPRIP_signaltonoise_kallisto_barplot.pdf',plot=p,device='pdf',width=12,height=5.5,unit='in') 

#########################################
pratio<-read.csv('Kot1_CLIPCLAPRIP_signaltonoise_kallisto_means.csv',header=T)

pratio  <- pratio  %>% mutate(
  condition = factor(condition, levels=c('RIP_mm10_HNRNPK','RIP_mm10_HNRNPU',
                                         'CLIP_hg38_HNRNPU','CLIP_mm10_HNRNPU',
                                         'CLAP_hg38_HNRNPU','CLAP_mm10_HNRNPU')),
  groups = factor(groups,levels=c('RBPlessCK','CK'))
)

pos_bar <- position_dodge(width = 0.8)

p<-ggplot() +
  geom_col(
    data = pratio,
    aes(
      x = condition,
      y = rpmsum,
      fill  = groups,                           
      group = interaction(condition, groups)       
    ),
    position = pos_bar,
    width    = 0.8,
    linewidth = 1.0
  ) +
  labs(x = 'Samples', y = "RPM Sum") +
  scale_fill_manual(values = c('#4f4a8c','#9d9d9d'))+
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


ggsave('Kot1_CLIPCLAPRIP_signaltonoise_kallisto_barplot.pdf',plot=p,device='pdf',width=12,height=5.5,unit='in') 

