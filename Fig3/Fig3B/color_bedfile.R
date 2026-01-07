# organize bedfiles

##########################################
# X
##########################################

filelist<-dir(pattern='Xist_.*_merged100.bed',full.names = F)


colordf<-data.frame(com=gsub('_neg_merged100.bed','',filelist),
                    color=c("240,100,142","70,183,73","54,196,234","94,110,105",
                            "176,153,200","190,156,47"))

if (exists('combbed')) {rm(combbed)}
for (f in filelist) {
  
  ft<-read.table(f, sep='\t',header=F)
  # extract community name
  cname<-gsub('_neg_merged100.bed','',f)
  # create 9 col bed with color
  ccol<-colordf$color[which(colordf$com==cname)]
  bname<-gsub('ist_','',cname)
  ft$V4<-bname
  
  ft$thickStart<-ft$V2
  ft$thickEnd<-ft$V3
  ft$itemRgb<-ccol
  
  if (exists('combbed')) {
    combbed<-rbind(combbed,ft)
  } else {
    combbed<-ft
  }
  
}

fnull<-read.table('Xist_null.bed', sep='\t',header=F)
fnull$V4<-'Xnull'
fnull$thickStart<-fnull$V2
fnull$thickEnd<-fnull$V3
fnull$itemRgb<-"236,236,236"

combbed<-rbind(combbed,fnull)

write.table(combbed,'Xist_commu_color.bed',row.names = F, col.names = F, quote = F)


##########################################
# A
##########################################

filelist<-dir(pattern='Airn_.*_merged100.bed',full.names = F)


colordf<-data.frame(com=gsub('_pos_merged100.bed','',filelist),
                    color=c("240,100,142","85,186,71","52,197,243","106,104,104",
                            "194,148,195","244,118,42","28,185,151","184,152,48"))

if (exists('combbed')) {rm(combbed)}
for (f in filelist) {
  
  ft<-read.table(f, sep='\t',header=F)
  # extract community name
  cname<-gsub('_pos_merged100.bed','',f)
  # create 9 col bed with color
  ccol<-colordf$color[which(colordf$com==cname)]
  bname<-gsub('irn_','',cname)
  ft$V4<-bname
  
  ft$thickStart<-ft$V2
  ft$thickEnd<-ft$V3
  ft$itemRgb<-ccol
  
  if (exists('combbed')) {
    combbed<-rbind(combbed,ft)
  } else {
    combbed<-ft
  }
  
}

fnull<-read.table('Airn_null.bed', sep='\t',header=F)
fnull$V4<-'Anull'
fnull$thickStart<-fnull$V2
fnull$thickEnd<-fnull$V3
fnull$itemRgb<-"236,236,236"

combbed<-rbind(combbed,fnull)

write.table(combbed,'Airn_commu_color.bed',row.names = F, col.names = F, quote = F)

##########################################
# K
##########################################
library(gtools)

filelist<-dir(pattern='Kot1_.*_merged100.bed',full.names = F)
filelist<-mixedsort(filelist)

colordf<-data.frame(com=gsub('_neg_merged100.bed','',filelist),
                    color=c("240,92,140","62,199,243","96,187,70","76,68,60",
                            "211,142,190","189,155,47","26,182,135","92,156,212",
                            "197,183,159"))

if (exists('combbed')) {rm(combbed)}
for (f in filelist) {
  
  ft<-read.table(f, sep='\t',header=F)
  # extract community name
  cname<-gsub('_neg_merged100.bed','',f)
  # create 9 col bed with color
  ccol<-colordf$color[which(colordf$com==cname)]
  bname<-gsub('ot1_','',cname)
  ft$V4<-bname
  
  ft$thickStart<-ft$V2
  ft$thickEnd<-ft$V3
  ft$itemRgb<-ccol
  
  if (exists('combbed')) {
    combbed<-rbind(combbed,ft)
  } else {
    combbed<-ft
  }
  
}

fnull<-read.table('Kot1_null.bed', sep='\t',header=F)
fnull$V4<-'Knull'
fnull$thickStart<-fnull$V2
fnull$thickEnd<-fnull$V3
fnull$itemRgb<-"236,236,236"

combbed<-rbind(combbed,fnull)

write.table(combbed,'Kot1_commu_color.bed',row.names = F, col.names = F, quote = F)

