# rank peaks based on CLAP data of meanr1r2rpmlessck*(meanr1r2rpmlessck/peakrpm)
# get top 10k, 5k 2.5k 1k and XIST and OT1 peaks
# for the selected peaks, calculate rpkm that is rpm normalized by peak length kb
# sum the rpkm for each selection group
# then use the selected peak bedfile to get counts in CLIP and calculate rpkm sum
# compare the summed rpkm of human vs mouse

htag_clap<-read.csv('SAFA_PlusTag_CLAP_hg38_allinfo.csv',header=T)
mtag_clap<-read.csv('SAFA_MinusTag_CLAP_mm10_allinfo.csv',header=T)

# calculate replicates 

htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm<-(htag_clap$SAFA_PlusTag_CLAP_Rep1_hg38_rpm+htag_clap$SAFA_PlusTag_CLAP_Rep2_hg38_rpm)/2
htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm_lessCK<-htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm-htag_clap$SAFA_PlusTag_CLAP_hg38_ck_rpm

htag_clap$rank_rpm<-htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm_lessCK*(htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm_lessCK/htag_clap$SAFA_PlusTag_CLAP_meanRep12_hg38_rpm)
mtag_clap$rank_rpm<-mtag_clap$SAFA_MinusTag_CLAP_mm10_rpm_lessCK*(mtag_clap$SAFA_MinusTag_CLAP_mm10_rpm_lessCK/mtag_clap$SAFA_MinusTag_CLAP_mm10_rpm)

# get the true peaks in the right organism and order based on rank_rpm

htag_clap_fl<-htag_clap[htag_clap$true_peaks,] # 72167
htag_clap_fl<-htag_clap_fl[!grepl('_mm10',htag_clap_fl$Chr),] # 72022
mtag_clap_fl<-mtag_clap[mtag_clap$true_peaks,] # 155122
mtag_clap_fl<-mtag_clap_fl[!grepl('_hg38',mtag_clap_fl$Chr),] #80883

htag_clap_fl<-htag_clap_fl[order(htag_clap_fl$rank_rpm,decreasing = T),]
mtag_clap_fl<-mtag_clap_fl[order(mtag_clap_fl$rank_rpm,decreasing = T),]

# calculate rpkm for each peak
htag_clap_fl$rpkm<-htag_clap_fl$SAFA_PlusTag_CLAP_hg38_rpm*1000/htag_clap_fl$Length
mtag_clap_fl$rpkm<-mtag_clap_fl$SAFA_MinusTag_CLAP_mm10_rpm*1000/mtag_clap_fl$Length
# and rpkm for rep1 and rep2
htag_clap_fl$SAFA_PlusTag_CLAP_Rep1_hg38_rpm_lessCK<-htag_clap_fl$SAFA_PlusTag_CLAP_Rep1_hg38_rpm-htag_clap_fl$SAFA_PlusTag_CLAP_hg38_ck_rpm
htag_clap_fl$SAFA_PlusTag_CLAP_Rep2_hg38_rpm_lessCK<-htag_clap_fl$SAFA_PlusTag_CLAP_Rep2_hg38_rpm-htag_clap_fl$SAFA_PlusTag_CLAP_hg38_ck_rpm

htag_clap_fl$rpkm_Rep1<-htag_clap_fl$SAFA_PlusTag_CLAP_Rep1_hg38_rpm*1000/htag_clap_fl$Length
htag_clap_fl$rpkm_Rep2<-htag_clap_fl$SAFA_PlusTag_CLAP_Rep2_hg38_rpm*1000/htag_clap_fl$Length

htag_clap_fl$rpkm_hpeak_mdata<-htag_clap_fl$SAFA_PlusTag_CLAP_hg38_ck_rpm*1000/htag_clap_fl$Length
mtag_clap_fl$rpkm_mpeak_hdata<-mtag_clap_fl$SAFA_MinusTag_CLAP_mm10_ck_rpm*1000/mtag_clap_fl$Length


htag_clap_fl$ranking<-seq(from=1, to=nrow(htag_clap_fl),by=1)
mtag_clap_fl$ranking<-seq(from=1, to=nrow(mtag_clap_fl),by=1)

write.csv(htag_clap_fl,'hg38CLAP_all_rpkm.csv',row.names = F)
write.csv(mtag_clap_fl,'mm10CLAP_all_rpkm.csv',row.names = F)


############## fliter and generate bedfile 
# human
# XIST chrX:73820656-73852714
# KOT1: chr11:2608328-2699994


(xidx<-which(htag_clap_fl$Chr=='chrX' & htag_clap_fl$Strand=='-' & htag_clap_fl$Start>=73820656 & htag_clap_fl$End<73852723))
(kidx<-which(htag_clap_fl$Chr=='chr11' & htag_clap_fl$Strand=='-' & htag_clap_fl$Start>=2608328 & htag_clap_fl$End<2699994))
combidx<-union(c(1:10000),c(xidx,kidx))
h_10xk<-htag_clap_fl[combidx,]


write.table(h_10xk[,c(28,2,3,4,5)],'hg38_CLAP_10xk.saf',
            row.names = F, col.names = F, quote = F,sep='\t')

############## initiate and fill data sheets
fstats<-matrix(nrow=72,ncol=4)
fstats<-as.data.frame(fstats)
colnames(fstats)<-c('condition','replicate','group','rpkmsum')
fstats$condition<-rep(c('hg38_CLAP','hg38_CLAP','hg38_CLAP','hg38_CLAP','hg38_CLIP','hg38_CLIP','hg38_CLIP','hg38_CLIP','mm10_CLAP','mm10_CLAP','mm10_CLIP','mm10_CLIP'),each=6)
fstats$replicate<-rep(c('Rep1','Rep2','mpeak_hdata_Rep1','mpeak_hdata_Rep2','Rep1','Rep2','mpeak_hdata_Rep1','mpeak_hdata_Rep2','Rep1','hpeak_mdata','Rep1','hpeak_mdata'),each=6)
fstats$group<-rep(c('top10k','top5k','top2.5k','top1k','XIST','KCNQ1OT1'),times=12)


(xidx<-which(htag_clap_fl$Chr=='chrX' & htag_clap_fl$Strand=='-' & htag_clap_fl$Start>=73820656 & htag_clap_fl$End<73852723))
(kidx<-which(htag_clap_fl$Chr=='chr11' & htag_clap_fl$Strand=='-' & htag_clap_fl$Start>=2608328 & htag_clap_fl$End<2699994))
h10k<-htag_clap_fl[c(1:10000),]
hx<-htag_clap_fl[xidx,]
hk<-htag_clap_fl[kidx,]

write.csv(htag_clap_fl[union(c(1:10000),c(xidx,kidx)),],'hg38CLAP_10xk_rpkm.csv',row.names = F)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_Rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_Rep1[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_Rep1[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_Rep1[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_Rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_Rep1)




fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_Rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_Rep2[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_Rep2[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_Rep2[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_Rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_Rep2)



fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_hpeak_mdata)

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_hpeak_mdata[1:5000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_hpeak_mdata[1:2500])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_hpeak_mdata[1:1000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_hpeak_mdata)

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_hpeak_mdata)


############### mouse
# fliter and generate bedfile and data sheets
# Xist chrX:103460366-103483254
# Kot1: chr7:143203458-143296549
(xidx<-which(mtag_clap_fl$Chr=='chrX' &mtag_clap_fl$Strand=='-' & mtag_clap_fl$Start>=103460366 & mtag_clap_fl$End<103483254))
(kidx<-which(mtag_clap_fl$Chr=='chr7' &mtag_clap_fl$Strand=='-' & mtag_clap_fl$Start>=143203458 & mtag_clap_fl$End<143296549))
combidx<-union(c(1:10000),c(xidx,kidx))
m_10xk<-mtag_clap_fl[combidx,]

write.table(m_10xk[,c(16,2,3,4,5)],'mm10_CLAP_10xk.saf',
            row.names = F, col.names = F, quote = F,sep='\t')

######################
# go to CLIP_fc.sh to count CLIP data as well as replicates data
#####################


######## load in replicates info
mmck1_10k<-read.table('mm10CLAP10k_CLAP_ck1_fc',sep='\t',header=T)
mmck2_10k<-read.table('mm10CLAP10k_CLAP_ck2_fc',sep='\t',header=T)

sum(mmck1_10k$Start!=m_10xk$Start)
sum(mmck2_10k$Start!=m_10xk$Start)


m_10xk$SAFA_MinusTag_CLAP_mm10_ck1_count<-mmck1_10k$SAFA_PlusTag_CLAP_Rep1_mm10_pairedonly.sam
m_10xk$SAFA_MinusTag_CLAP_mm10_ck2_count<-mmck2_10k$SAFA_PlusTag_CLAP_Rep2_mm10_pairedonly.sam


m_10xk$SAFA_MinusTag_CLAP_mm10_ck1_rpm<-m_10xk$SAFA_MinusTag_CLAP_mm10_ck1_count*1000000/23065288
m_10xk$SAFA_MinusTag_CLAP_mm10_ck2_rpm<-m_10xk$SAFA_MinusTag_CLAP_mm10_ck2_count*1000000/21401736


m_10xk$rpkm_mpeak_hdata_rep1<-m_10xk$SAFA_MinusTag_CLAP_mm10_ck1_rpm*1000/m_10xk$Length
m_10xk$rpkm_mpeak_hdata_rep2<-m_10xk$SAFA_MinusTag_CLAP_mm10_ck2_rpm*1000/m_10xk$Length

write.csv(m_10xk,'mm10CLAP_10xk_rpkm.csv',row.names = F)

(xidx<-which(m_10xk$Chr=='chrX' &m_10xk$Strand=='-' & m_10xk$Start>=103460366 & m_10xk$End<103483254))
(kidx<-which(m_10xk$Chr=='chr7' &m_10xk$Strand=='-' & m_10xk$Start>=143203458 & m_10xk$End<143296549))

m10k<-m_10xk[c(1:10000),]
mx<-m_10xk[xidx,]
mk<-m_10xk[kidx,]


fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm)

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm[1:5000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm[1:2500])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm[1:1000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='XIST')]<-sum(mx$rpkm)

fstats$rpkmsum[which(fstats$condition=='mm10_CLAP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm)




fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm_mpeak_hdata_rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='XIST')]<-sum(mx$rpkm_mpeak_hdata_rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm_mpeak_hdata_rep1)




fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm_mpeak_hdata_rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='XIST')]<-sum(mx$rpkm_mpeak_hdata_rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLAP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm_mpeak_hdata_rep2)


################### CLIP using CLAP defined peak bedfile
# in CLIP_fc.sh to use featurecounts to count CLIP under CLAP defined peaks

htag_clip_10k<-read.table('hg38CLAP10k_CLIP_fc',sep='\t',header=T)
htag_clip_10k_ck<-read.table('hg38CLAP10k_CLIP_ck_fc',sep='\t',header=T)
htag_clip_10k_Rep1<-read.table('hg38CLAP10k_CLIP_Rep1_fc',sep='\t',header=T)
htag_clip_10k_Rep2<-read.table('hg38CLAP10k_CLIP_Rep2_fc',sep='\t',header=T)
sum(htag_clip_10k$Start!=htag_clip_10k_ck$Start)
sum(htag_clip_10k$End!=htag_clip_10k_Rep1$End)
sum(htag_clip_10k$End!=htag_clip_10k_Rep2$End)


colnames(htag_clip_10k)[7]<-'hg38CLAP10k_CLIP_count'
htag_clip_10k$hg38CLAP10k_CLIP_Rep1_count<-htag_clip_10k_Rep1$SAFA_PlusTag_CLIP_Rep1_hg38_pairedonly.sam
htag_clip_10k$hg38CLAP10k_CLIP_Rep2_count<-htag_clip_10k_Rep2$SAFA_PlusTag_CLIP_Rep2_hg38_pairedonly.sam
htag_clip_10k$hg38CLAP10k_CLIP_ck_count<-htag_clip_10k_ck$SAFA_MinusTag_CLIP_hg38_pairedonly.sam

# calculate rpm by counts*1000000/totalfastqreads and rpkm
htag_clip_10k$hg38CLAP10k_CLIP_rpm<-htag_clip_10k$hg38CLAP10k_CLIP_count*1000000/38511502
htag_clip_10k$hg38CLAP10k_CLIP_Rep1_rpm<-htag_clip_10k$hg38CLAP10k_CLIP_Rep1_count*1000000/20105697
htag_clip_10k$hg38CLAP10k_CLIP_Rep2_rpm<-htag_clip_10k$hg38CLAP10k_CLIP_Rep2_count*1000000/18405805
htag_clip_10k$hg38CLAP10k_CLIP_ck_rpm<-htag_clip_10k$hg38CLAP10k_CLIP_ck_count*1000000/4879320

htag_clip_10k$hg38CLAP10k_CLIP_rpm_lessCK<-htag_clip_10k$hg38CLAP10k_CLIP_rpm-htag_clip_10k$hg38CLAP10k_CLIP_ck_rpm
htag_clip_10k$hg38CLAP10k_CLIP_Rep1_rpm_lessCK<-htag_clip_10k$hg38CLAP10k_CLIP_Rep1_rpm-htag_clip_10k$hg38CLAP10k_CLIP_ck_rpm
htag_clip_10k$hg38CLAP10k_CLIP_Rep2_rpm_lessCK<-htag_clip_10k$hg38CLAP10k_CLIP_Rep2_rpm-htag_clip_10k$hg38CLAP10k_CLIP_ck_rpm


htag_clip_10k$rpkm<-htag_clip_10k$hg38CLAP10k_CLIP_rpm*1000/htag_clip_10k$Length
htag_clip_10k$rpkm_Rep1<-htag_clip_10k$hg38CLAP10k_CLIP_Rep1_rpm*1000/htag_clip_10k$Length
htag_clip_10k$rpkm_Rep2<-htag_clip_10k$hg38CLAP10k_CLIP_Rep2_rpm*1000/htag_clip_10k$Length
htag_clip_10k$rpkm_hpeak_mdata<-htag_clip_10k$hg38CLAP10k_CLIP_ck_rpm*1000/htag_clip_10k$Length


write.csv(htag_clip_10k,'hg38CLAP10xk_CLIP_rpkm.csv',row.names = F)

(xidx<-which(htag_clip_10k$Chr=='chrX' & htag_clip_10k$Strand=='-' & htag_clip_10k$Start>=73820656 & htag_clip_10k$End<73852723))
(kidx<-which(htag_clip_10k$Chr=='chr11' & htag_clip_10k$Strand=='-' & htag_clip_10k$Start>=2608328 & htag_clip_10k$End<2699994))
sum(htag_clip_10k$Start!=h_10xk$Start)

h10k<-htag_clip_10k[c(1:10000),]
hx<-htag_clip_10k[xidx,]
hk<-htag_clip_10k[kidx,]


fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_Rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_Rep1[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_Rep1[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_Rep1[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_Rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_Rep1)




fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_Rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_Rep2[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_Rep2[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_Rep2[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_Rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='Rep2' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_Rep2)




fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top10k')]<-sum(h10k$rpkm_hpeak_mdata)

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top5k')]<-sum(h10k$rpkm_hpeak_mdata[1:5000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top2.5k')]<-sum(h10k$rpkm_hpeak_mdata[1:2500])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='top1k')]<-sum(h10k$rpkm_hpeak_mdata[1:1000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='XIST')]<-sum(hx$rpkm_hpeak_mdata)

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='hpeak_mdata' &
                       fstats$group=='KCNQ1OT1')]<-sum(hk$rpkm_hpeak_mdata)

######### mouse
mtag_clip_10k<-read.table('mm10CLAP10k_CLIP_fc',sep='\t',header=T)
mtag_clip_10k_ck<-read.table('mm10CLAP10k_CLIP_ck_fc',sep='\t',header=T)
mtag_clip_10k_ck1<-read.table('mm10CLAP10k_CLIP_ck1_fc',sep='\t',header=T)
mtag_clip_10k_ck2<-read.table('mm10CLAP10k_CLIP_ck2_fc',sep='\t',header=T)
sum(mtag_clip_10k$Start!=mtag_clip_10k_ck$Start)
sum(mtag_clip_10k$Start!=mtag_clip_10k_ck1$Start)
sum(mtag_clip_10k$Start!=mtag_clip_10k_ck2$Start)



colnames(mtag_clip_10k)[7]<-'mm10CLAP10k_CLIP_count'
mtag_clip_10k$mm10CLAP10k_CLIP_ck_count<-mtag_clip_10k_ck$SAFA_PlusTag_CLIP_mm10_pairedonly.sam
mtag_clip_10k$mm10CLAP10k_CLIP_ck1_count<-mtag_clip_10k_ck1$SAFA_PlusTag_CLIP_Rep1_mm10_pairedonly.sam
mtag_clip_10k$mm10CLAP10k_CLIP_ck2_count<-mtag_clip_10k_ck2$SAFA_PlusTag_CLIP_Rep2_mm10_pairedonly.sam


# calculate rpm by counts*1000000/totalfastqreads and rpkm
mtag_clip_10k$mm10CLAP10k_CLIP_rpm<-mtag_clip_10k$mm10CLAP10k_CLIP_count*1000000/4879320
mtag_clip_10k$mm10CLAP10k_CLIP_ck_rpm<-mtag_clip_10k$mm10CLAP10k_CLIP_ck_count*1000000/38511502
mtag_clip_10k$mm10CLAP10k_CLIP_rpm_lessCK<-mtag_clip_10k$mm10CLAP10k_CLIP_rpm-mtag_clip_10k$mm10CLAP10k_CLIP_ck_rpm
mtag_clip_10k$mm10CLAP10k_CLIP_ck1_rpm<-mtag_clip_10k$mm10CLAP10k_CLIP_ck1_count*1000000/20105697
mtag_clip_10k$mm10CLAP10k_CLIP_ck2_rpm<-mtag_clip_10k$mm10CLAP10k_CLIP_ck2_count*1000000/18405805


mtag_clip_10k$rpkm<-mtag_clip_10k$mm10CLAP10k_CLIP_rpm*1000/mtag_clip_10k$Length
mtag_clip_10k$rpkm_mpeak_hdata<-mtag_clip_10k$mm10CLAP10k_CLIP_ck_rpm*1000/mtag_clip_10k$Length
mtag_clip_10k$rpkm_mpeak_hdata_rep1<-mtag_clip_10k$mm10CLAP10k_CLIP_ck1_rpm*1000/mtag_clip_10k$Length
mtag_clip_10k$rpkm_mpeak_hdata_rep2<-mtag_clip_10k$mm10CLAP10k_CLIP_ck2_rpm*1000/mtag_clip_10k$Length


write.csv(mtag_clip_10k,'mm10CLAP10xk_CLIP_rpkm.csv',row.names = F)


(xidx<-which(mtag_clip_10k$Chr=='chrX' &mtag_clip_10k$Strand=='-' & mtag_clip_10k$Start>=103460366 & mtag_clip_10k$End<103483254))
(kidx<-which(mtag_clip_10k$Chr=='chr7' &mtag_clip_10k$Strand=='-' & mtag_clip_10k$Start>=143203458 & mtag_clip_10k$End<143296549))
sum(mtag_clip_10k$Start!=m_10xk$Start)
m10k<-mtag_clip_10k[c(1:10000),]
mx<-mtag_clip_10k[xidx,]
mk<-mtag_clip_10k[kidx,]



fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm)

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm[1:5000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm[1:2500])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm[1:1000])

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='XIST')]<-sum(mx$rpkm)

fstats$rpkmsum[which(fstats$condition=='mm10_CLIP' & 
                       fstats$replicate=='Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm)





fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm_mpeak_hdata_rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm_mpeak_hdata_rep1[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='XIST')]<-sum(mx$rpkm_mpeak_hdata_rep1)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep1' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm_mpeak_hdata_rep1)





fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top10k')]<-sum(m10k$rpkm_mpeak_hdata_rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top5k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:5000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top2.5k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:2500])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='top1k')]<-sum(m10k$rpkm_mpeak_hdata_rep2[1:1000])

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='XIST')]<-sum(mx$rpkm_mpeak_hdata_rep2)

fstats$rpkmsum[which(fstats$condition=='hg38_CLIP' & 
                       fstats$replicate=='mpeak_hdata_Rep2' &
                       fstats$group=='KCNQ1OT1')]<-sum(mk$rpkm_mpeak_hdata_rep2)



write.csv(fstats,'CLIP_CLAP_rpkm_stats.csv',row.names = F)
