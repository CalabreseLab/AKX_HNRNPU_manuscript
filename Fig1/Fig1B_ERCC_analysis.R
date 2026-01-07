
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
## e.g., counts <- as.matrix(read.csv("raw_counts.csv", row.names=1))

combc<-comb[,c(7:71)]
rownames(combc)<-comb$GeneID
# 65 samples

# 2) drop features with all zeros
combc <- combc[rowSums(combc > 0) > 0, , drop=FALSE]
# 65 samples

# 3) compute the geometric mean per feature across samples,
#    ignoring zeros (DESeq2’s behavior)
gm <- apply(combc, 1, function(v) {
  vpos <- v[v > 0]
  if (length(vpos) == 0) return(0)
  exp(mean(log(vpos)))
})

# keep only features with positive geometric mean
# keep <- gm > 0
# counts2 <- counts[keep, , drop=FALSE]
# gm <- gm[keep]

sum(gm>0) # 92

# 4) for each sample, compute ratios feature_count / feature_geomean
# row-wise margin=1
ratios <- sweep(combc, 1, gm, "/")

# 5) sample size factor = median of those ratios (ignoring NAs and zeros)
sf <- apply(ratios, 2, median, na.rm = TRUE)

sf<-as.data.frame(sf)
prop<-rownames(sf)
prop[which(prop=='nxf1_counts')]<-'nxf1_rip1_counts'
prop[which(prop=='u2af65_counts')]<-'u2af65_rip1_counts'

prop<-gsub('_counts','',prop)
prop<-strsplit(prop,'_')

sf$RBP<-sapply(prop,'[[',1)
sf$replicate<-sapply(prop,'[[',2)
colnames(sf)[1]<-'mor'
# ensure sample is a factor (keeps desired order on x-axis)
sf$RBP<-toupper(sf$RBP)
sf$RBP[which(sf$RBP=='IGG')]<-'IgG'

new_levels <- c('IgG', setdiff(unique(sf$RBP), 'IgG'))

sf$RBP <- factor(sf$RBP, levels = new_levels)

# Compute per-sample means
sf_mean <- sf %>%
  group_by(RBP) %>%
  summarise(mean_mor = mean(mor), .groups = "drop")

# normalized MoR
# plot for each replicate: [avg mor of IgG]/[mor of replicate]

sf_mean$IgGnorm_mean_mor<-sf_mean$mean_mor[which(sf_mean$RBP=='IgG')]/sf_mean$mean_mor

sf$IgGnorm_mor<-sf_mean$mean_mor[which(sf_mean$RBP=='IgG')]/sf$mor

p<-ggplot() +
  geom_hline(yintercept = 1,color='#a50000',linewidth=1.0)+
  geom_col(data = sf_mean,aes(x = RBP, y = IgGnorm_mean_mor),
           width = 0.6,fill = "grey80", color='grey80') +
  geom_point(data = sf,aes(x = RBP, y = IgGnorm_mor),
             position = position_jitter(width = 0.12, height = 0),
             size = 2.2,alpha = 0.8) +
  scale_y_log10()+
  labs(x = NULL, y = "IgG normed MoR") +
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major=element_line(color='grey',linewidth=0.1),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.title=element_blank(),
        legend.position='none',
        axis.text.x=element_text(size=26,color='black',angle=45,hjust=1),
        axis.text.y=element_text(size=26,color='black'),
        axis.title.x=element_text(size=26),
        axis.title.y=element_text(size=26,angle=90))

ggsave('ERCC_Meidan_of_ratios_barplot_IgGnorm.pdf',plot=p,device='pdf',width=13,height=5,unit='in') 


