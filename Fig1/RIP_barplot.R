# plot RIP HNRNPK and HNRNPU reassco ratios

library(ggplot2)
library(dplyr)

pratio<-read.csv('RIP_plot_ratios_Fig1G.csv',header=T)
ppoint<-read.csv('RIP_plot_points_Fig1G.csv',header=T)

pratio  <- pratio  %>% mutate(
  file = factor(file, levels=c('GFP','HNRNPK','HNRNPU')),
  exps  = factor(exps,levels=c('paired_reads','peaks')),
  methods = factor(methods,levels=c('Input','RIP','top10k','top5k','top2.5k','top1k','XIST','KCNQ1OT1'))
)
ppoint  <- ppoint  %>% mutate(
  file = factor(file, levels=levels(pratio$file)),
  exps  = factor(exps,  levels = levels(pratio$exps)),
  methods = factor(methods, levels = levels(pratio$methods))
)

## use a common dodging position so bars and points align
pos_dodge <- position_dodge2(width = 0.8,preserve = 'single',padding = 0)

p<-ggplot() +
  # BAR LAYER: one bar per row in pratio
  geom_col(
    data = pratio,
    aes(
      x = file,
      y = normalized_perct,
      fill  = methods,                           
      color = exps,                          
      group = interaction(exps, methods)       # defines dodging groups
    ),
    position = pos_dodge,
    width    = 0.8,
    linewidth = 1.0
  ) +
  
  # POINT LAYER: replicate points from ppoint over the same bars
  geom_point(
    data = ppoint,
    aes(
      x = file,
      y = normalized_perct,                          
      group = interaction(exps, methods)       # same grouping for same dodge
    ),
    color = 'black',
    position = pos_dodge,                      # EXACTLY the same position object
    size     = 3.0,
    alpha    = 0.6
  ) +
  # some typical cosmetics
  labs(x = 'Samples', y = "Percentage") +
  scale_fill_manual(values = c('#1b9e77','#d95f02','#393564','#4f4a8c','#655faa','#8581bc','#a5a2ce','#c5c3e0'))+
  scale_color_manual(values = c('#1e1e1e','#d8d8d8'))+
  theme(panel.background=element_rect(fill='transparent'),
        plot.margin = margin(t=85, r=5, b=5, l=5, "pt"),
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


ggsave('RIP_Fig1G.pdf',plot=p,device='pdf',width=13,height=4.35,unit='in') 

