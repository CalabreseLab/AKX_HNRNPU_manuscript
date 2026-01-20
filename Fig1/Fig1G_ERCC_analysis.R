setwd("/Users/shuang/CalabreseLab/RIP_data_process/hkhu_wiggle/ERCC")

library(ggplot2)
library(dplyr)


# Compute a per-sample size factor from ERCCs raw counts using the median-of-ratios
# Intuition: for each ERCC, take the ratio of its count in a sample to the geometric mean across samples
# the median of those ratios is the size factor. 
# This uses all ERCCs, downweights outliers, and doesn’t require absolute calibration.

filelist<-dir(pattern='*._fc_rpm.txt')

for (f in filelist) {
  
  df<-read.table(f,header=T,sep='\t')
  
  counts<-df[,grepl('counts',colnames(df)),drop=F]
  
  if (exists('comb')) {
    if(sum(df$GeneID!=comb$GeneID)==0){
      comb<-cbind(comb,counts)
    } else {
      print('matching error')
      print(f)
      break
    }
  } else {
    comb<-cbind(df[,c(1:6)],counts)
  }
  
}

## counts: integer matrix/data.frame of raw counts
## rows = features (ERCCs if using spike-ins), cols = samples

combc<-comb[,c(7:30)]
rownames(combc)<-comb$GeneID
# 6 samples

# drop features with all zeros
combc <- combc[rowSums(combc > 0) > 0, , drop=FALSE]
# 90 ERCC

# compute the geometric mean per feature (ERCC) across samples,
# ignoring zeros (DESeq2’s behavior)
gm <- apply(combc, 1, function(v) {
  vpos <- v[v > 0]
  if (length(vpos) == 0) return(0)
  exp(mean(log(vpos)))
})

# keep only features with positive geometric mean
# keep <- gm > 0
# counts2 <- counts[keep, , drop=FALSE]
# gm <- gm[keep]

sum(gm>0) # 90

# 4) for each sample, compute ratios feature_count / feature_geomean
# row-wise margin=1
ratios <- sweep(combc, 1, gm, "/")

# 5) sample size factor = median of those ratios (ignoring NAs and zeros)
sf <- apply(ratios, 2, median, na.rm = TRUE)

sf<-as.data.frame(sf)
prop<-rownames(sf)

prop<-gsub('_counts','',prop)
prop<-strsplit(prop,'_rip')

sf$RBP<-sapply(prop,'[[',1)
sf$replicate<-sapply(prop,'[[',2)
colnames(sf)[1]<-'mor'
# ensure sample is a factor (keeps desired order on x-axis)

new_levels<-c('hm_GFP_flag','hm_HNRNPK_flag','hm_HNRNPU_flag')

sf$RBP <- factor(sf$RBP, levels = new_levels)

# Compute per-sample means
sf_mean <- sf %>%
  group_by(RBP) %>%
  summarise(mean_mor = mean(mor), .groups = "drop")

#######################################
# plot for each replicate: [avg mor of IgG]/[mor of replicate]

sf_mean$norm_mean_mor<-NA
sf_mean$norm_mean_mor<-sf_mean$mean_mor[which(sf_mean$RBP=='hm_GFP_flag')]/sf_mean$mean_mor


sf$norm_mor<-NA
sf$norm_mor<-sf_mean$mean_mor[which(sf_mean$RBP=='hm_GFP_flag')]/sf$mor

p<-ggplot() +
  geom_hline(yintercept = 1,color='#a50000',linewidth=1.0)+
  geom_col(data = sf_mean,aes(x = RBP, y = norm_mean_mor),
           width = 0.6,fill = "grey80", color='grey80') +
  geom_point(data = sf,aes(x = RBP, y = norm_mor),
             position = position_jitter(width = 0.12, height = 0),
             size = 2.2,alpha = 0.8) +
  scale_y_log10()+
  labs(x = NULL, y = "Normalized\nMoR") +
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=25, r=5, b=5, l=5, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='none',
        #legend.key.height = unit(1,'line'),
        #legend.key.width  = unit(1,'line'),
        #legend.key=element_rect(fill='transparent'),
        #legend.text=element_text(size=26),
        #legend.text.align=0,
        axis.text.x=element_text(size=26,color='black',angle=45,hjust=1),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))

ggsave('ERCC_Meidan_of_ratios_barplot_norm.pdf',plot=p,device='pdf',width=3,height=4.2,unit='in') 



