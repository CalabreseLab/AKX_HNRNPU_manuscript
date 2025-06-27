# dealing with all edges including negative ones
x<-read.csv('Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3_edges.csv',header=T)
a<-read.csv('Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced_edges.csv',header=T)
k<-read.csv('Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced_edges.csv',header=T)

#####################
# calculate inter and intra community r values for XAK, plot violin plot and calculate p values
x$r <- ifelse(x$Color == "Negative", -x$Weight, x$Weight)
a$r <- ifelse(a$Color == "Negative", -a$Weight, a$Weight)
k$r <- ifelse(k$Color == "Negative", -k$Weight, k$Weight)

xnode<-read.csv('Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3_nodes_ld.csv',header=T)
anode<-read.csv('Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced_nodes_ld.csv',header=T)
knode<-read.csv('Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced_nodes_ld.csv',header=T)

add_commu<-function(df, dfnode) {
  
  df$community<-''
  
  for (n in 1:nrow(df)) {
    ts<-x$Source[n]
    tt<-x$Target[n]
    if (dfnode$Color[which(dfnode$Id==ts)]==dfnode$Color[which(dfnode$Id==tt)]) {
      df$community[n]<-'intra'
    } else {
      df$community[n]<-'inter'
    }
  }
  
  
  return(df)
  
}

x<-add_commu(x,xnode)
a<-add_commu(a,anode)
k<-add_commu(k,knode)

library(ggplot2)

# data def not normal dist
# use wilcox to test for sig
wilcox.test(x$r[which(x$community=='intra')],x$r[which(x$community=='inter')],alternative='greater')
wilcox.test(a$r[which(a$community=='intra')],a$r[which(a$community=='inter')],alternative='greater')
wilcox.test(k$r[which(k$community=='intra')],k$r[which(k$community=='inter')],alternative='greater')
# all pvals: p-value < 2.2e-16

# compare XAK intra
wilcox.test(x$r[which(x$community=='intra')],a$r[which(a$community=='intra')],alternative='less')
# 0.3573
wilcox.test(x$r[which(x$community=='intra')],k$r[which(k$community=='intra')],alternative='less')
# 0.03951
wilcox.test(a$r[which(a$community=='intra')],k$r[which(k$community=='intra')],alternative='less')
# 0.006577

# compare XAK inter
wilcox.test(x$r[which(x$community=='inter')],a$r[which(a$community=='inter')],alternative='less')
# 5.359e-12
wilcox.test(x$r[which(x$community=='inter')],k$r[which(k$community=='inter')],alternative='less')
# < 2.2e-16
wilcox.test(a$r[which(a$community=='inter')],k$r[which(k$community=='inter')],alternative='less')
# < 2.2e-16


# violin plot to show the shape better
# combine XAK data together

x$gene<-'Xist'
a$gene<-'Airn'
k$gene<-'Kcnq1ot1'

comb<-rbind(x,a,k)
comb$gene <- factor(comb$gene, levels = c("Xist", "Airn", "Kcnq1ot1"))

p<-ggplot(data=comb,aes(x = community, y = r)) +
  geom_violin(aes(fill=community,color=community)) +
  labs(x = "Community",y = "Pearson's r Value") +
  theme(plot.title=element_text(size=22),
        panel.background=element_rect(fill='white'),
        panel.grid.major=element_line(color='grey',linewidth =0.3),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        legend.position='none',
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size=18,angle=45, hjust=1),
        axis.text.y=element_text(size=18),
        strip.text = element_text(size = 20)) +
  facet_wrap(~gene,ncol=3)+
  scale_fill_manual(values=c('#af8dc3','#7fbf7b'))+
  scale_color_manual(values=c('#762a83','#1b7837'))

ggsave("XAK_inter_intra_violin.pdf", plot = p, width = 4.5, height = 4, units = "in")
