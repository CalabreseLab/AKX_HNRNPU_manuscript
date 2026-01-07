###############################################
# Make histograms and csv table of all transcripts RIP signal Pearson r values relative to XAK
# draw vertical lines marking XAK
###############################################
library(ggplot2)
library(ggrepel)
xist<-"Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3"
airn<-"Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced"
ot1<-"Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"

newdf<-read.csv('TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024_filtered_07222024.csv',header=T)

newdf_dat<-newdf[,c(7:33)]
row.names(newdf_dat)<-newdf$target_id
xist_dat<-newdf_dat[xist,]
airn_dat<-newdf_dat[airn,]
ot1_dat<-newdf_dat[ot1,]

xist_rval <- apply(newdf_dat, 1, function(row) {
  cor(as.numeric(row), as.numeric(xist_dat), method = "pearson")
})


airn_rval <- apply(newdf_dat, 1, function(row) {
  cor(as.numeric(row), as.numeric(airn_dat), method = "pearson")
})

ot1_rval <- apply(newdf_dat, 1, function(row) {
  cor(as.numeric(row), as.numeric(ot1_dat), method = "pearson")
})

comb_rval<-as.data.frame(newdf$target_id)
colnames(comb_rval)<-'target_id'
comb_rval$xist_rval<-as.numeric(xist_rval)
comb_rval$airn_rval<-as.numeric(airn_rval)
comb_rval$ot1_rval<-as.numeric(ot1_rval)

write.csv(comb_rval,'XAK_pearson_r_RIP.csv',row.names = F)

comb_rval$fillcol<-'A'

###########################################
# create plots with XAK_pearson_r_RIP.csv

#################xist

x_lnc<-c(airn,ot1)
x_lnc_names<-c('Airn','Kcnq1ot1')

x_dat<-comb_rval$xist_rval[match(x_lnc,comb_rval$target_id)]

x_lncdat<-as.data.frame(x_lnc_names)
colnames(x_lncdat)<-'label'
x_lncdat$rval<-x_dat

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = comb_rval, aes(x = xist_rval, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=4000, by=1000)) +
  labs(x = "All vs Xist R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.4,0.8), ylim=c(0,4000))+
  scale_x_continuous(breaks=seq(from=-0.4, to=0.8, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=24))

p<-p + geom_vline(data = x_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = x_lncdat, aes(x = rval, y = 3800, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("allvsXist_hist_pearson_RIP.pdf", plot = p, width = 5.5, height = 5, units = "in")

#################airn

a_lnc<-c(xist,ot1)
a_lnc_names<-c('Xist','Kcnq1ot1')

a_dat<-comb_rval$airn_rval[match(a_lnc,comb_rval$target_id)]

a_lncdat<-as.data.frame(a_lnc_names)
colnames(a_lncdat)<-'label'
a_lncdat$rval<-a_dat

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = comb_rval, aes(x = airn_rval, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=2000, by=1000)) +
  labs(x = "All vs Airn R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.6,0.9), ylim=c(0,2000))+
  scale_x_continuous(breaks=seq(from=-0.6, to=0.9, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=24))

p<-p + geom_vline(data = a_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = a_lncdat, aes(x = rval, y = 1800, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -1,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=5)


ggsave("allvsAirn_hist_pearson_RIP.pdf", plot = p, width = 5.5, height = 5, units = "in")


#################ot1

k_lnc<-c(xist,airn)
k_lnc_names<-c('Xist','Airn')

k_dat<-comb_rval$ot1_rval[match(k_lnc,comb_rval$target_id)]

k_lncdat<-as.data.frame(k_lnc_names)
colnames(k_lncdat)<-'label'
k_lncdat$rval<-k_dat

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = comb_rval, aes(x = ot1_rval, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=2000, by=1000)) +
  labs(x = "All vs Kcnq1ot1 R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.6,0.9), ylim=c(0,2000))+
  scale_x_continuous(breaks=seq(from=-0.6, to=0.9, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x=element_text(size=24, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=24))

p<-p + geom_vline(data = k_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = k_lncdat, aes(x = rval, y = 1800, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = 0, vjust = -1,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps = 20,
                  force=5)


ggsave("allvsKcnq1ot1_hist_pearson_RIP.pdf", plot = p, width = 5.5, height = 5, units = "in")
