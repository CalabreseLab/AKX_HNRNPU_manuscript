# quantify 27 protein RBP network similarities between all transcripts

#################################################
# Make histograms and csv table of all 19295 transcripts network Pearson r values relative to XAK
# replace NAs in the edge wrights with 0
#################################################
library(dplyr)
library(data.table)

dt<-fread('comb_allgenes_cor_08092024_long.csv',header=T)
# 6772545/19295=351
# remove genes that contain NAs in the weight for better correlation calculation
# genes_with_na <- dt[is.na(Weight), unique(gene)]
# 5464 genes
# Remove rows that belong to genes with NA in the weight column
# dt_clean <- dt[!gene %in% genes_with_na]
# 4854681 rows
# unique genes 13831

dt$Weight[is.na(dt$Weight)]<-0


xist<-'Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3'
airn<-'Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced'
ot1<-'Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced'


# loop thru all genes against XAK


xak <- data.table(
  gene = character(),
  xist_r = numeric(),
  airn_r = numeric(),
  kcnq1ot1_r = numeric()
)

genelist<-unique(dt$gene)

# Filter rows for each gene using data.table's syntax
xdt <- dt[gene == xist]
adt <- dt[gene == airn]
kdt <- dt[gene == ot1]

# check if xak RBP pairs are in the same order
all(xdt$RBP1 == adt$RBP1 & xdt$RBP2 == adt$RBP2)
all(xdt$RBP1 == kdt$RBP1 & xdt$RBP2 == kdt$RBP2)

for (n in 1:length(genelist)) {
  
  if (n %% 100 ==0) {print(n)}
  
  g<-genelist[n]
  
  dt2 <- dt[gene == g]
  
  # Check if the RBP1 and RBP2 are the same
  if (all(xdt$RBP1 == dt2$RBP1 & xdt$RBP2 == dt2$RBP2)) {
    xcor <- cor(xdt$Weight, dt2$Weight, method = "pearson")
    acor <- cor(adt$Weight, dt2$Weight, method = "pearson")
    kcor <- cor(kdt$Weight, dt2$Weight, method = "pearson")
  } else {
    print('RBP1 and RBP2 matching problem')
    break
  }
  
  
  
  temp <- data.table(
    gene = g,
    xist_r = xcor,
    airn_r = acor,
    kcnq1ot1_r = kcor
  )
  
  xak<-rbindlist(list(xak, temp))
  
}

# still NA exists, as there are genes with all NAs as the edge weights before
# right now is all 0, and cor will procude NA with 0 sd
# check NA num: 72, unique genes 19223
sum(is.na(xak$xist_r)) # 72
sum(!is.na(xak$airn_r)) # 19223
sum(!is.na(xak$kcnq1ot1_r))
# list all NA genes
xak$gene[is.na(xak$xist_r)]


write.csv(xak,'XAK_vs_all_pearson_r_network_NA0.csv',row.names = F)

###############
# plot hist for all vs X/A/K with XAK + others labeled
library(ggplot2)
library(ggrepel)


xak$fillcol<-'A'

#### xist ####

x_lnc<-c(airn,ot1)
x_lnc_names<-c('Airn','Kcnq1ot1')

x_dat<-xak$xist_r[match(x_lnc,xak$gene)]

x_lncdat<-as.data.frame(x_lnc_names)
colnames(x_lncdat)<-'label'
x_lncdat$rval<-x_dat

# remove NAs from the table
x_lncdat<-x_lncdat[!is.na(x_lncdat$rval),]

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = xak, aes(x = xist_r, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=3500, by=1000)) +
  labs(x = "All vs Xist R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.2,0.75), ylim=c(0,3500))+
  scale_x_continuous(breaks=seq(from=-0.2, to=0.75, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.x=element_text(size=26, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=26))

p<-p + geom_vline(data = x_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = x_lncdat, aes(x = rval, y = 3350, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("allvsXist_hist_pearson_network_NA0.pdf", plot = p, width = 5.5, height = 5, units = "in")


#### airn ####

a_lnc<-c(xist,ot1)
a_lnc_names<-c('Xist','Kcnq1ot1')

a_dat<-xak$airn_r[match(a_lnc,xak$gene)]

a_lncdat<-as.data.frame(a_lnc_names)
colnames(a_lncdat)<-'label'
a_lncdat$rval<-a_dat

# remove NAs from the table
a_lncdat<-a_lncdat[!is.na(a_lncdat$rval),]

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = xak, aes(x = airn_r, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=3500, by=1000)) +
  labs(x = "All vs Airn R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.2,0.75), ylim=c(0,3500))+
  scale_x_continuous(breaks=seq(from=-0.2, to=0.75, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.x=element_text(size=26, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=26))

p<-p + geom_vline(data = a_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = a_lncdat, aes(x = rval, y = 3350, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("allvsAirn_hist_pearson_network_NA0.pdf", plot = p, width = 5.5, height = 5, units = "in")


#### ot1 ####

k_lnc<-c(xist,airn)
k_lnc_names<-c('Xist','Airn')

k_dat<-xak$kcnq1ot1_r[match(k_lnc,xak$gene)]

k_lncdat<-as.data.frame(k_lnc_names)
colnames(k_lncdat)<-'label'
k_lncdat$rval<-k_dat

# remove NAs from the table
k_lncdat<-k_lncdat[!is.na(k_lncdat$rval),]

hbreaks<-seq(from=-1, to=1, by=0.05)
# Create the plot
p<-ggplot() +
  geom_histogram(data = xak, aes(x = kcnq1ot1_r, y = after_stat(count), fill=fillcol,color=fillcol), breaks=hbreaks, alpha=0.5) +
  scale_y_continuous(name = "Count",breaks=seq(from=0, to=3500, by=1000)) +
  labs(x = "All vs Kcnq1ot1 R Values",y = "Count") +
  coord_cartesian(xlim=c(-0.2,0.75), ylim=c(0,3500))+
  scale_x_continuous(breaks=seq(from=-0.2, to=0.75, by=0.2),labels = scales::number_format(accuracy = 0.1))+
  scale_color_manual(values=c('#5a5a5a','#2c7bb6'))+
  scale_fill_manual(values=c('grey','#abd9e9'))+
  theme(panel.background=element_rect(fill='white'),
        panel.grid.major=element_blank(),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.x=element_text(size=26, angle = 45, vjust=0.5),
        axis.text.y=element_text(size=26))

p<-p + geom_vline(data = k_lncdat, aes(xintercept = rval), color = "#1b9e77", linetype = "solid",linewidth=0.5)+
  geom_text_repel(data = k_lncdat, aes(x = rval, y = 3350, label = label), 
                  size = 6, nudge_x = 0, direction = "y", hjust = -0.05, vjust = -0.3,
                  segment.color = "darkred",   # Color of the connecting line
                  segment.size = 0.5,       # Thickness of the connecting line
                  segment.linetype = "solid",
                  max.overlaps=20,
                  force=3)


ggsave("allvsKcnq1ot1_hist_pearson_network_NA0.pdf", plot = p, width = 5.5, height = 5, units = "in")

