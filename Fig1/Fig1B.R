# plot heatmaps of rankings and top1k pearsons

library(ggplot2)
library(ComplexHeatmap)
library(circlize)


########### ranking heatmap #################
df<-read.csv('TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024.csv',header=T)
# 114783 in total
# filter filter for length >=500, median exp >0.0625, chrom-fraction >0.75, 
df<-df[which(df$length.x >= 500 & df$Median_tpm >0.0625 & df$chrom_enrichment >0.75),]
# 19295

# Retain length.x,Median_tpm,chrom_average,cyto_average,chrom_enrichment
# Retain all rips except bmi1,ezh2,suz12,epop,mtf2,jarid2
# retain “X_over_igg” column (and rename X_less_igg)
newdf<-df[,c(1,2,22,34,35,36,seq(from=38,to=104,by=2))]
newdf<-newdf[,c(1:7,9,12:16,19:21,23:34,36:40)]

colnames(newdf)<-gsub('rpm_over_igg','rpm_less_igg',colnames(newdf),fixed=T)

write.csv(newdf,'TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024_filtered_07222024.csv',row.names = F)


rankdf<-as.data.frame(matrix(nrow=29,ncol=3))
colnames(rankdf)<-c('Xist','Airn','Kcnq1ot1')
rownames(rankdf)<-c('expression_total','expression_chr',gsub('_rpm_less_igg','',colnames(newdf)[7:33]))

######expression total median tpm
df_sorted <- newdf[order(newdf[,'Median_tpm'], decreasing = TRUE), ]

# Get the rank of specific genes 
rankdf[1,]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                    "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                    "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)

######expression chr
df_sorted <- newdf[order(newdf[,'chrom_average'], decreasing = TRUE), ]

# Get the rank of specific genes 
rankdf[2,]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                    "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                    "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)


#### rbps
rbplist<-colnames(newdf)[7:33]

for (n in 1:length(rbplist)) {
  
  tcol<-rbplist[n]
  df_sorted <- newdf[order(newdf[,tcol], decreasing = TRUE), ]
  
  rankdf[(2+n),]<-match(c("Xist_chrX_103460366_103483254_-_ENSMUSG00000086503.4_ENSMUST00000127786.3", 
                          "Airn_chr17_12741398_12830151_+_ENSMUSG00000078247.5_ENSMUSG00000078247.5.unspliced", 
                          "Kcnq1ot1_chr7_143203458_143296549_-_ENSMUSG00000101609.2_ENSMUST00000185789.2.monoexonic.unspliced"),df_sorted$target_id)
  
  
}

########### plot

# change the order to AKX
rankdf<-rankdf[,c('Airn','Kcnq1ot1','Xist')]

log2_rankdf <- log2(rankdf)
# Create a color mapping from white to purple
color_mapping <- colorRamp2(c(min(log2_rankdf), max(log2_rankdf)), c("#34009a", "white"))

cn = colnames(rankdf)

# Create labels for the legend
legend_labels <- c('1','8','64','512','4096',as.character(nrow(newdf)))
legend_values <- c(1,8,64,512,4096,nrow(newdf))
# Transform these values to match the scale of the heatmap
transformed_legend_values <- log2(legend_values)


log2_rankdf<-as.matrix(log2_rankdf)
row.names(log2_rankdf)<-toupper(row.names(log2_rankdf))


# Open a PDF device
pdf("ranking_heatmap_01112025.pdf", width = 5.5, height = 12)  # Adjust 'width' and 'height' as needed

# Plot the heatmap
ht<-Heatmap(log2_rankdf, 
            col = color_mapping, 
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            rect_gp = gpar(col = "#e0e0e0", lwd = 1),
            show_row_names = TRUE,
            show_column_names = FALSE,
            #column_names_side = "top",
            row_names_side = "left",
            #column_names_gp = gpar(fontsize = 18,fontface = 3),
            row_names_gp = gpar(fontsize = 18),
            top_annotation = HeatmapAnnotation(
              text = anno_text(cn, rot = 0, location = unit(0.1, "npc"), just = "center",
                               gp = gpar(fontsize = 18, fontface = "italic")),
              annotation_height = max_text_width(cn, gp = gpar(fontsize = 18,fontface = "italic"))
            ),
            cell_fun = function(j, i, x, y, width, height, fill) {
              text_color <- ifelse(rankdf[i, j] <= 150, "white", "black")
              grid.text(rankdf[i, j], x, y, gp = gpar(col = text_color, fontsize = 18))
            },
            heatmap_legend_param = list(
              title = "Ranking among\nchromatin enriched transcripts",
              at = transformed_legend_values,
              labels = legend_labels,
              color_bar = "continuous",
              legend_direction = "horizontal",
              legend_width = unit(8, "cm"),
              #legend_height = unit(2, "cm"),
              labels_gp = gpar(fontsize = 18),
              title_gp = gpar(fontsize = 18))
        
)

draw(ht, heatmap_legend_side = "top",padding = unit(c(2, 5, 2, 2), "mm"))

dev.off()

