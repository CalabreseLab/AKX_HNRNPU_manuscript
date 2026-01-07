# organzie featurecounts results and locate true peaks
# > 2x ck in 2 replicates or in the only replicate

############### human peaks##############################

############ HK
fc<-read.table('h_HNRNPK_merg_fc',sep='\t',header=T)
rep1<-read.table('h_HNRNPK_r1_fc',sep='\t',header=T)
rep2<-read.table('h_HNRNPK_r2_fc',sep='\t',header=T)
ck<-read.table('h_HNRNPK_ck_fc',sep='\t',header=T)

sum(fc$Start!=rep1$Start)
sum(fc$End!=rep2$End)
sum(fc$Length!=ck$Length)

colnames(fc)[7]<-'h_HNRNPK_merg_count'

fc$h_HNRNPK_r1_count<-rep1$h_HNRNPK_r1_Aligned.out.sam
fc$h_HNRNPK_r2_count<-rep2$h_HNRNPK_r2_Aligned.out.sam
fc$h_HNRNPK_ck_count<-ck$h_IGG_Aligned.out.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$h_HNRNPK_merg_rpm<-fc$h_HNRNPK_merg_count*1000000/54031188
fc$h_HNRNPK_r1_rpm<-fc$h_HNRNPK_r1_count*1000000/28438612
fc$h_HNRNPK_r2_rpm<-fc$h_HNRNPK_r2_count*1000000/25592576
fc$h_HNRNPK_ck_rpm<-fc$h_HNRNPK_ck_count*1000000/35715333
fc$h_HNRNPK_r12mean_rpm<-(fc$h_HNRNPK_r1_rpm+fc$h_HNRNPK_r2_rpm)/2

fc$h_HNRNPK_r12mean_lessCK_rpm<-fc$h_HNRNPK_r12mean_rpm-fc$h_HNRNPK_ck_rpm

# check if each peak is > 2x ck

fc$rep1_true<-(fc$h_HNRNPK_r1_rpm>2*fc$h_HNRNPK_ck_rpm)
fc$rep2_true<-(fc$h_HNRNPK_r2_rpm>2*fc$h_HNRNPK_ck_rpm)

fc$true_peaks<-(fc$rep1_true & fc$rep2_true)

write.csv(fc,'h_HNRNPK_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
#fctrue<-fctrue[!grepl('_mm10',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$h_HNRNPK_r12mean_lessCK_rpm*(fctrue$h_HNRNPK_r12mean_lessCK_rpm/fctrue$h_HNRNPK_r12mean_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,22,20,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'h_HNRNPK_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')


############ HU
fc<-read.table('h_HNRNPU_merg_fc',sep='\t',header=T)
rep1<-read.table('h_HNRNPU_r1_fc',sep='\t',header=T)
rep2<-read.table('h_HNRNPU_r2_fc',sep='\t',header=T)
ck<-read.table('h_HNRNPU_ck_fc',sep='\t',header=T)

sum(fc$Start!=rep1$Start)
sum(fc$End!=rep2$End)
sum(fc$Length!=ck$Length)

colnames(fc)[7]<-'h_HNRNPU_merg_count'

fc$h_HNRNPU_r1_count<-rep1$h_HNRNPU_r1_Aligned.out.sam
fc$h_HNRNPU_r2_count<-rep2$h_HNRNPU_r2_Aligned.out.sam
fc$h_HNRNPU_ck_count<-ck$h_IGG_Aligned.out.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$h_HNRNPU_merg_rpm<-fc$h_HNRNPU_merg_count*1000000/69648708
fc$h_HNRNPU_r1_rpm<-fc$h_HNRNPU_r1_count*1000000/34022770
fc$h_HNRNPU_r2_rpm<-fc$h_HNRNPU_r2_count*1000000/35625938
fc$h_HNRNPU_ck_rpm<-fc$h_HNRNPU_ck_count*1000000/35715333
fc$h_HNRNPU_r12mean_rpm<-(fc$h_HNRNPU_r1_rpm+fc$h_HNRNPU_r2_rpm)/2

fc$h_HNRNPU_r12mean_lessCK_rpm<-fc$h_HNRNPU_r12mean_rpm-fc$h_HNRNPU_ck_rpm

# check if each peak is > 2x ck

fc$rep1_true<-(fc$h_HNRNPU_r1_rpm>2*fc$h_HNRNPU_ck_rpm)
fc$rep2_true<-(fc$h_HNRNPU_r2_rpm>2*fc$h_HNRNPU_ck_rpm)

fc$true_peaks<-(fc$rep1_true & fc$rep2_true)

write.csv(fc,'h_HNRNPU_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
#fctrue<-fctrue[!grepl('_mm10',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$h_HNRNPU_r12mean_lessCK_rpm*(fctrue$h_HNRNPU_r12mean_lessCK_rpm/fctrue$h_HNRNPU_r12mean_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,22,20,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'h_HNRNPU_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')

############### mouse hm peaks##############################

############ HK
fc<-read.table('hm_HNRNPK_flag_merg_fc',sep='\t',header=T)
rep1<-read.table('hm_HNRNPK_flag_r1_fc',sep='\t',header=T)
rep2<-read.table('hm_HNRNPK_flag_r2_fc',sep='\t',header=T)
ck<-read.table('hm_HNRNPK_flag_ck_fc',sep='\t',header=T)

sum(fc$Start!=rep1$Start)
sum(fc$End!=rep2$End)
sum(fc$Length!=ck$Length)

colnames(fc)[7]<-'hm_HNRNPK_flag_merg_count'

fc$hm_HNRNPK_flag_r1_count<-rep1$hm_HNRNPK_flag_r1_Aligned.out.sam
fc$hm_HNRNPK_flag_r2_count<-rep2$hm_HNRNPK_flag_r2_Aligned.out.sam
fc$hm_HNRNPK_flag_ck_count<-ck$hm_GFP_flag_Aligned.out.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$hm_HNRNPK_flag_merg_rpm<-fc$hm_HNRNPK_flag_merg_count*1000000/62978874
fc$hm_HNRNPK_flag_r1_rpm<-fc$hm_HNRNPK_flag_r1_count*1000000/30746853
fc$hm_HNRNPK_flag_r2_rpm<-fc$hm_HNRNPK_flag_r2_count*1000000/32232021
fc$hm_HNRNPK_flag_ck_rpm<-fc$hm_HNRNPK_flag_ck_count*1000000/60319397
fc$hm_HNRNPK_flag_r12mean_rpm<-(fc$hm_HNRNPK_flag_r1_rpm+fc$hm_HNRNPK_flag_r2_rpm)/2

fc$hm_HNRNPK_flag_r12mean_lessCK_rpm<-fc$hm_HNRNPK_flag_r12mean_rpm-fc$hm_HNRNPK_flag_ck_rpm

# check if each peak is > 2x ck

fc$rep1_true<-(fc$hm_HNRNPK_flag_r1_rpm>2*fc$hm_HNRNPK_flag_ck_rpm)
fc$rep2_true<-(fc$hm_HNRNPK_flag_r2_rpm>2*fc$hm_HNRNPK_flag_ck_rpm)

fc$true_peaks<-(fc$rep1_true & fc$rep2_true)

write.csv(fc,'hm_HNRNPK_flag_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
fctrue<-fctrue[!grepl('_hg38',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$hm_HNRNPK_flag_r12mean_lessCK_rpm*(fctrue$hm_HNRNPK_flag_r12mean_lessCK_rpm/fctrue$hm_HNRNPK_flag_r12mean_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,22,20,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'hm_HNRNPK_flag_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')


############ HU
fc<-read.table('hm_HNRNPU_flag_merg_fc',sep='\t',header=T)
rep1<-read.table('hm_HNRNPU_flag_r1_fc',sep='\t',header=T)
rep2<-read.table('hm_HNRNPU_flag_r2_fc',sep='\t',header=T)
ck<-read.table('hm_HNRNPU_flag_ck_fc',sep='\t',header=T)

sum(fc$Start!=rep1$Start)
sum(fc$End!=rep2$End)
sum(fc$Length!=ck$Length)

colnames(fc)[7]<-'hm_HNRNPU_flag_merg_count'

fc$hm_HNRNPU_flag_r1_count<-rep1$hm_HNRNPU_flag_r1_Aligned.out.sam
fc$hm_HNRNPU_flag_r2_count<-rep2$hm_HNRNPU_flag_r2_Aligned.out.sam
fc$hm_HNRNPU_flag_ck_count<-ck$hm_GFP_flag_Aligned.out.sam

# calculate rpm by counts*1000000/totalfastqreads
# totalfastqreads comes from preivous counting results
fc$hm_HNRNPU_flag_merg_rpm<-fc$hm_HNRNPU_flag_merg_count*1000000/69381631
fc$hm_HNRNPU_flag_r1_rpm<-fc$hm_HNRNPU_flag_r1_count*1000000/32401416
fc$hm_HNRNPU_flag_r2_rpm<-fc$hm_HNRNPU_flag_r2_count*1000000/36980215
fc$hm_HNRNPU_flag_ck_rpm<-fc$hm_HNRNPU_flag_ck_count*1000000/60319397
fc$hm_HNRNPU_flag_r12mean_rpm<-(fc$hm_HNRNPU_flag_r1_rpm+fc$hm_HNRNPU_flag_r2_rpm)/2

fc$hm_HNRNPU_flag_r12mean_lessCK_rpm<-fc$hm_HNRNPU_flag_r12mean_rpm-fc$hm_HNRNPU_flag_ck_rpm

# check if each peak is > 2x ck

fc$rep1_true<-(fc$hm_HNRNPU_flag_r1_rpm>2*fc$hm_HNRNPU_flag_ck_rpm)
fc$rep2_true<-(fc$hm_HNRNPU_flag_r2_rpm>2*fc$hm_HNRNPU_flag_ck_rpm)

fc$true_peaks<-(fc$rep1_true & fc$rep2_true)

write.csv(fc,'hm_HNRNPU_flag_allinfo.csv',row.names = F)

fctrue<-fc[fc$true_peaks,]
fctrue<-fctrue[!grepl('_hg38',fctrue$Chr),]
fctrue$score<-0

fctrue$rank_rpm<-fctrue$hm_HNRNPU_flag_r12mean_lessCK_rpm*(fctrue$hm_HNRNPU_flag_r12mean_lessCK_rpm/fctrue$hm_HNRNPU_flag_r12mean_rpm)
fctrue<-fctrue[order(fctrue$rank_rpm,decreasing = T),]
fctrue$ranking<-seq(from=1, to=nrow(fctrue),by=1)
fctrue$color<-'197,197,197'
fctrue$color[which(fctrue$ranking<=1000)]<-'128,0,128'
fctrue$color[which(fctrue$ranking<=5000 & fctrue$ranking>1000)]<-'246,0,0'
fctrue$color[which(fctrue$ranking<=10000 & fctrue$ranking>5000)]<-'209,69,255'
fctrue$color[which(fctrue$ranking<=20000 & fctrue$ranking>10000)]<-'255,192,203'


# write bedfile and color the regions based on ranking

fcbed<-fctrue[,c(2,3,4,22,20,5)]
fcbed$Start<-fcbed$Start-1
fcbed$thickstart<-fctrue$Start-1
fcbed$thickend<-fctrue$End
fcbed$color<-fctrue$color

write.table(fcbed,'hm_HNRNPU_flag_truepeaks.bed',
            row.names = F, col.names = F, quote = F,sep='\t')


