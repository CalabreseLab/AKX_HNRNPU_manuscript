# organzie featurecounts results and locate true peaks
# > 2x ck in 2 replicates or in the only replicate
# featurecounts -p count the read pairs (fragments)

#############################################
fc<-read.table('SAFA_PlusTag_CLAP_hg38_fc',sep='\t',header=T)
rep1<-read.table('SAFA_PlusTag_CLAP_Rep1_hg38_fc',sep='\t',header=T)
rep2<-read.table('SAFA_PlusTag_CLAP_Rep2_hg38_fc',sep='\t',header=T)
ck<-read.table('SAFA_PlusTag_CLAP_hg38_ck_fc',sep='\t',header=T)

sum(fc$Start!=rep1$Start)
sum(fc$End!=rep2$End)
sum(fc$Length!=ck$Length)

colnames(fc)[7]<-'SAFA_PlusTag_CLAP_hg38_count'

fc$SAFA_PlusTag_CLAP_Rep1_hg38_count<-rep1$SAFA_PlusTag_CLAP_Rep1_hg38_pairedonly.sam
fc$SAFA_PlusTag_CLAP_Rep2_hg38_count<-rep2$SAFA_PlusTag_CLAP_Rep2_hg38_pairedonly.sam
fc$SAFA_PlusTag_CLAP_hg38_ck_count<-ck$SAFA_MinusTag_CLAP_hg38_pairedonly.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$SAFA_PlusTag_CLAP_hg38_rpm<-fc$SAFA_PlusTag_CLAP_hg38_count*1000000/44467024
fc$SAFA_PlusTag_CLAP_Rep1_hg38_rpm<-fc$SAFA_PlusTag_CLAP_Rep1_hg38_count*1000000/23065288
fc$SAFA_PlusTag_CLAP_Rep2_hg38_rpm<-fc$SAFA_PlusTag_CLAP_Rep2_hg38_count*1000000/21401736
fc$SAFA_PlusTag_CLAP_hg38_ck_rpm<-fc$SAFA_PlusTag_CLAP_hg38_ck_count*1000000/5703926
fc$SAFA_PlusTag_CLAP_hg38_rpm_lessCK<-fc$SAFA_PlusTag_CLAP_hg38_rpm-fc$SAFA_PlusTag_CLAP_hg38_ck_rpm

# check if each peak is > 2x ck

fc$rep1_true<-(fc$SAFA_PlusTag_CLAP_Rep1_hg38_rpm>2*fc$SAFA_PlusTag_CLAP_hg38_ck_rpm)
fc$rep2_true<-(fc$SAFA_PlusTag_CLAP_Rep2_hg38_rpm>2*fc$SAFA_PlusTag_CLAP_hg38_ck_rpm)

fc$true_peaks<-(fc$rep1_true & fc$rep2_true)

write.csv(fc,'SAFA_PlusTag_CLAP_hg38_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
fctrue<-fctrue[!grepl('_mm10',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$SAFA_PlusTag_CLAP_hg38_rpm_lessCK*(fctrue$SAFA_PlusTag_CLAP_hg38_rpm_lessCK/fctrue$SAFA_PlusTag_CLAP_hg38_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,21,19,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'SAFA_PlusTag_CLAP_hg38_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')


############################################

fc<-read.table('SAFA_MinusTag_CLAP_mm10_fc',sep='\t',header=T)
ck<-read.table('SAFA_MinusTag_CLAP_mm10_ck_fc',sep='\t',header=T)

sum(fc$Start!=ck$Start)

colnames(fc)[7]<-'SAFA_MinusTag_CLAP_mm10_count'

fc$SAFA_MinusTag_CLAP_mm10_ck_count<-ck$SAFA_PlusTag_CLAP_mm10_pairedonly.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$SAFA_MinusTag_CLAP_mm10_rpm<-fc$SAFA_MinusTag_CLAP_mm10_count*1000000/5703926
fc$SAFA_MinusTag_CLAP_mm10_ck_rpm<-fc$SAFA_MinusTag_CLAP_mm10_ck_count*1000000/44467024
fc$SAFA_MinusTag_CLAP_mm10_rpm_lessCK<-fc$SAFA_MinusTag_CLAP_mm10_rpm-fc$SAFA_MinusTag_CLAP_mm10_ck_rpm

# check if each peak is > 2x ck

fc$true_peaks<-(fc$SAFA_MinusTag_CLAP_mm10_rpm>2*fc$SAFA_MinusTag_CLAP_mm10_ck_rpm)

write.csv(fc,'SAFA_MinusTag_CLAP_mm10_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
fctrue<-fctrue[!grepl('_hg38',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$SAFA_MinusTag_CLAP_mm10_rpm_lessCK*(fctrue$SAFA_MinusTag_CLAP_mm10_rpm_lessCK/fctrue$SAFA_MinusTag_CLAP_mm10_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,15,13,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'SAFA_MinusTag_CLAP_mm10_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')


