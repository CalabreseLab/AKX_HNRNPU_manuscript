# rank peaks based on meanr1r2rpmlessck*(meanr1r2rpmlessck/meanr1r2rpm)
# get top 10k, 5k 2.5k 1k and XIST peaks in human and in mouse
# in mouse for the selected peaks, calculate rpkm that is rpm (not normalized) by peak length kb
# sum the rpkm for each selection group
# then use the human peaks to get human counts in hm samples
# compare the summed rpkm of human vs mouse

m_hk<-read.csv('hm_HNRNPK_flag_allinfo.csv',header=T)
m_hu<-read.csv('hm_HNRNPU_flag_allinfo.csv',header=T)

h_hk<-read.csv('h_HNRNPK_allinfo.csv',header=T)
h_hu<-read.csv('h_HNRNPU_allinfo.csv',header=T)

# calculate rank rpm
m_hk$rank_rpm<-m_hk$hm_HNRNPK_flag_r12mean_lessCK_rpm*(m_hk$hm_HNRNPK_flag_r12mean_lessCK_rpm/m_hk$hm_HNRNPK_flag_r12mean_rpm)
m_hu$rank_rpm<-m_hu$hm_HNRNPU_flag_r12mean_lessCK_rpm*(m_hu$hm_HNRNPU_flag_r12mean_lessCK_rpm/m_hu$hm_HNRNPU_flag_r12mean_rpm)

h_hk$rank_rpm<-h_hk$h_HNRNPK_r12mean_lessCK_rpm*(h_hk$h_HNRNPK_r12mean_lessCK_rpm/h_hk$h_HNRNPK_r12mean_rpm)
h_hu$rank_rpm<-h_hu$h_HNRNPU_r12mean_lessCK_rpm*(h_hu$h_HNRNPU_r12mean_lessCK_rpm/h_hu$h_HNRNPU_r12mean_rpm)


# get the true peaks in the right organism and order based on rank_rpm

m_hk_fl<-m_hk[m_hk$true_peaks,]
m_hk_fl<-m_hk_fl[!grepl('_hg38',m_hk_fl$Chr),]

m_hu_fl<-m_hu[m_hu$true_peaks,]
m_hu_fl<-m_hu_fl[!grepl('_hg38',m_hu_fl$Chr),]

h_hk_fl<-h_hk[h_hk$true_peaks,]
h_hu_fl<-h_hu[h_hu$true_peaks,]

m_hk_fl<-m_hk_fl[order(m_hk_fl$rank_rpm,decreasing = T),]
m_hu_fl<-m_hu_fl[order(m_hu_fl$rank_rpm,decreasing = T),]
h_hk_fl<-h_hk_fl[order(h_hk_fl$rank_rpm,decreasing = T),]
h_hu_fl<-h_hu_fl[order(h_hu_fl$rank_rpm,decreasing = T),]


# and rpkm for rep1 and rep2, no normalized by ck - used for bar plot
# rpm normalize by ck - used for signal to noise
m_hk_fl$hm_HNRNPK_flag_r1_rpm_lessCK<-m_hk_fl$hm_HNRNPK_flag_r1_rpm-m_hk_fl$hm_HNRNPK_flag_ck_rpm
m_hk_fl$hm_HNRNPK_flag_r2_rpm_lessCK<-m_hk_fl$hm_HNRNPK_flag_r2_rpm-m_hk_fl$hm_HNRNPK_flag_ck_rpm

sum(m_hk_fl$hm_HNRNPK_flag_r1_rpm_lessCK<0)
sum(m_hk_fl$hm_HNRNPK_flag_r2_rpm_lessCK<0)

m_hk_fl$rpkm_r1<-m_hk_fl$hm_HNRNPK_flag_r1_rpm*1000/m_hk_fl$Length
m_hk_fl$rpkm_r2<-m_hk_fl$hm_HNRNPK_flag_r2_rpm*1000/m_hk_fl$Length


m_hu_fl$hm_HNRNPU_flag_r1_rpm_lessCK<-m_hu_fl$hm_HNRNPU_flag_r1_rpm-m_hu_fl$hm_HNRNPU_flag_ck_rpm
m_hu_fl$hm_HNRNPU_flag_r2_rpm_lessCK<-m_hu_fl$hm_HNRNPU_flag_r2_rpm-m_hu_fl$hm_HNRNPU_flag_ck_rpm

sum(m_hu_fl$hm_HNRNPU_flag_r1_rpm_lessCK<0)
sum(m_hu_fl$hm_HNRNPU_flag_r2_rpm_lessCK<0)


m_hu_fl$rpkm_r1<-m_hu_fl$hm_HNRNPU_flag_r1_rpm*1000/m_hu_fl$Length
m_hu_fl$rpkm_r2<-m_hu_fl$hm_HNRNPU_flag_r2_rpm*1000/m_hu_fl$Length


h_hk_fl$h_HNRNPK_r1_rpm_lessCK<-h_hk_fl$h_HNRNPK_r1_rpm-h_hk_fl$h_HNRNPK_ck_rpm
h_hk_fl$h_HNRNPK_r2_rpm_lessCK<-h_hk_fl$h_HNRNPK_r2_rpm-h_hk_fl$h_HNRNPK_ck_rpm

sum(h_hk_fl$h_HNRNPK_r1_rpm_lessCK<0)
sum(h_hk_fl$h_HNRNPK_r2_rpm_lessCK<0)

h_hk_fl$rpkm_r1<-h_hk_fl$h_HNRNPK_r1_rpm*1000/h_hk_fl$Length
h_hk_fl$rpkm_r2<-h_hk_fl$h_HNRNPK_r2_rpm*1000/h_hk_fl$Length


h_hu_fl$h_HNRNPU_r1_rpm_lessCK<-h_hu_fl$h_HNRNPU_r1_rpm-h_hu_fl$h_HNRNPU_ck_rpm
h_hu_fl$h_HNRNPU_r2_rpm_lessCK<-h_hu_fl$h_HNRNPU_r2_rpm-h_hu_fl$h_HNRNPU_ck_rpm

sum(h_hu_fl$h_HNRNPU_r1_rpm_lessCK<0)
sum(h_hu_fl$h_HNRNPU_r2_rpm_lessCK<0)

h_hu_fl$rpkm_r1<-h_hu_fl$h_HNRNPU_r1_rpm*1000/h_hu_fl$Length
h_hu_fl$rpkm_r2<-h_hu_fl$h_HNRNPU_r2_rpm*1000/h_hu_fl$Length


m_hk_fl$ranking<-seq(from=1, to=nrow(m_hk_fl),by=1)
m_hu_fl$ranking<-seq(from=1, to=nrow(m_hu_fl),by=1)
h_hk_fl$ranking<-seq(from=1, to=nrow(h_hk_fl),by=1)
h_hu_fl$ranking<-seq(from=1, to=nrow(h_hu_fl),by=1)


write.csv(m_hk_fl,'mm10_HNRNPK_all_rpkm.csv',row.names = F)
write.csv(m_hu_fl,'mm10_HNRNPU_all_rpkm.csv',row.names = F)
write.csv(h_hk_fl,'hg38_HNRNPK_all_rpkm.csv',row.names = F)
write.csv(h_hu_fl,'hg38_HNRNPU_all_rpkm.csv',row.names = F)


############## fliter and generate bedfile 

# XIST chrX:73820656-73852723
# KOT1: chr11:2608328-2699994

(xidx<-which(h_hk_fl$Chr=='chrX' & h_hk_fl$Strand=='-' & h_hk_fl$Start>=73820656 & h_hk_fl$End<73852723))
(kidx<-which(h_hk_fl$Chr=='chr11' & h_hk_fl$Strand=='-' & h_hk_fl$Start>=2608328 & h_hk_fl$End<2699994))
combidx<-union(c(1:10000),c(xidx,kidx))
h_hk_10xk<-h_hk_fl[combidx,]
write.table(h_hk_10xk[,c(25,2,3,4,5)],'hg38_HNRNPK_10xk.saf',
            row.names = F, col.names = F, quote = F,sep='\t')

write.csv(h_hk_fl[union(c(1:10000),c(xidx,kidx)),],'hg38_HNRNPK_10xk_rpkm.csv',row.names = F)

(xidx<-which(h_hu_fl$Chr=='chrX' & h_hu_fl$Strand=='-' & h_hu_fl$Start>=73820656 & h_hu_fl$End<73852723))
(kidx<-which(h_hu_fl$Chr=='chr11' & h_hu_fl$Strand=='-' & h_hu_fl$Start>=2608328 & h_hu_fl$End<2699994))
combidx<-union(c(1:10000),c(xidx,kidx))
h_hu_10xk<-h_hu_fl[combidx,]
write.table(h_hu_10xk[,c(25,2,3,4,5)],'hg38_HNRNPU_10xk.saf',
            row.names = F, col.names = F, quote = F,sep='\t')

write.csv(h_hu_fl[union(c(1:10000),c(xidx,kidx)),],'hg38_HNRNPU_10xk_rpkm.csv',row.names = F)

#############################
# go to RIP_fc.sh to count reads for human peaks in hm_ samples
##############################

############### fill in data
fstats<-matrix(nrow=48,ncol=4)
fstats<-as.data.frame(fstats)
colnames(fstats)<-c('condition','replicate','group','rpkmsum')
fstats$condition<-rep(c('HNRNPK','HNRNPU'),each=24)
fstats$replicate<-rep(c('mouse_r1','mouse_r2','human_r1','human_r2','mouse_r1','mouse_r2','human_r1','human_r2'),each=6)
fstats$group<-rep(c('top10k','top5k','top2.5k','top1k','XIST','KCNQ1OT1'),times=8)


############################### mouse
# Xist: chrX:103460366-103483254
# Kot1: chr7:143203458-143296549


##### HK
(xidx<-which(m_hk_fl$Chr=='chrX' & m_hk_fl$Strand=='-' & m_hk_fl$Start>=103460366 & m_hk_fl$End<103483254))
(kidx<-which(m_hk_fl$Chr=='chr7' & m_hk_fl$Strand=='-' & m_hk_fl$Start>=143203458 & m_hk_fl$End<143296549))
m_hk_10k<-m_hk_fl[c(1:10000),]
m_hk_x<-m_hk_fl[xidx,]
m_hk_k<-m_hk_fl[kidx,]

write.csv(m_hk_fl[union(c(1:10000),c(xidx,kidx)),],'mm10_HNRNPK_10xk_rpkm.csv',row.names = F)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top10k')]<-sum(m_hk_10k$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top5k')]<-sum(m_hk_10k$rpkm_r1[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top2.5k')]<-sum(m_hk_10k$rpkm_r1[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top1k')]<-sum(m_hk_10k$rpkm_r1[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='XIST')]<-sum(m_hk_x$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='KCNQ1OT1')]<-sum(m_hk_k$rpkm_r1)



fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top10k')]<-sum(m_hk_10k$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top5k')]<-sum(m_hk_10k$rpkm_r2[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top2.5k')]<-sum(m_hk_10k$rpkm_r2[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top1k')]<-sum(m_hk_10k$rpkm_r2[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='XIST')]<-sum(m_hk_x$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='KCNQ1OT1')]<-sum(m_hk_k$rpkm_r2)



##### HU
(xidx<-which(m_hu_fl$Chr=='chrX' & m_hu_fl$Strand=='-' & m_hu_fl$Start>=103460366 & m_hu_fl$End<103483254))
(kidx<-which(m_hu_fl$Chr=='chr7' & m_hu_fl$Strand=='-' & m_hu_fl$Start>=143203458 & m_hu_fl$End<143296549))
m_hu_10k<-m_hu_fl[c(1:10000),]
m_hu_x<-m_hu_fl[xidx,]
m_hu_k<-m_hu_fl[kidx,]

write.csv(m_hu_fl[union(c(1:10000),c(xidx,kidx)),],'mm10_HNRNPU_10xk_rpkm.csv',row.names = F)


fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top10k')]<-sum(m_hu_10k$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top5k')]<-sum(m_hu_10k$rpkm_r1[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top2.5k')]<-sum(m_hu_10k$rpkm_r1[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='top1k')]<-sum(m_hu_10k$rpkm_r1[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='XIST')]<-sum(m_hu_x$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r1' &
                       fstats$group=='KCNQ1OT1')]<-sum(m_hu_k$rpkm_r1)



fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top10k')]<-sum(m_hu_10k$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top5k')]<-sum(m_hu_10k$rpkm_r2[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top2.5k')]<-sum(m_hu_10k$rpkm_r2[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='top1k')]<-sum(m_hu_10k$rpkm_r2[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='XIST')]<-sum(m_hu_x$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='mouse_r2' &
                       fstats$group=='KCNQ1OT1')]<-sum(m_hu_k$rpkm_r2)




############### human
# count hm samples with human defined peaks from h_ samples
# results from RIP_fc.sh


hmh_hk<-read.table('hm_HNRNPK_flag_merg_hg38_fc',sep='\t',header=T)
hmh_hk_r1<-read.table('hm_HNRNPK_flag_r1_hg38_fc',sep='\t',header=T)
hmh_hk_r2<-read.table('hm_HNRNPK_flag_r2_hg38_fc',sep='\t',header=T)
hmh_hk_ck<-read.table('hm_HNRNPK_flag_ck_hg38_fc',sep='\t',header=T)
sum(hmh_hk$Start!=hmh_hk_ck$Start)
sum(hmh_hk$End!=hmh_hk_r1$End)
sum(hmh_hk$End!=hmh_hk_r2$End)


hmh_hu<-read.table('hm_HNRNPU_flag_merg_hg38_fc',sep='\t',header=T)
hmh_hu_r1<-read.table('hm_HNRNPU_flag_r1_hg38_fc',sep='\t',header=T)
hmh_hu_r2<-read.table('hm_HNRNPU_flag_r2_hg38_fc',sep='\t',header=T)
hmh_hu_ck<-read.table('hm_HNRNPU_flag_ck_hg38_fc',sep='\t',header=T)
sum(hmh_hu$Start!=hmh_hu_ck$Start)
sum(hmh_hu$End!=hmh_hu_r1$End)
sum(hmh_hu$End!=hmh_hu_r2$End)


colnames(hmh_hk)[7]<-'hm_HNRNPK_flag_hg38_merg_count'
hmh_hk$hm_HNRNPK_flag_hg38_r1_count<-hmh_hk_r1$hm_HNRNPK_flag_r1_hg38_Aligned.out.sam
hmh_hk$hm_HNRNPK_flag_hg38_r2_count<-hmh_hk_r2$hm_HNRNPK_flag_r2_hg38_Aligned.out.sam
hmh_hk$hm_HNRNPK_flag_hg38_ck_count<-hmh_hk_ck$hm_GFP_flag_hg38_Aligned.out.sam


colnames(hmh_hu)[7]<-'hm_HNRNPU_flag_hg38_merg_count'
hmh_hu$hm_HNRNPU_flag_hg38_r1_count<-hmh_hu_r1$hm_HNRNPU_flag_r1_hg38_Aligned.out.sam
hmh_hu$hm_HNRNPU_flag_hg38_r2_count<-hmh_hu_r2$hm_HNRNPU_flag_r2_hg38_Aligned.out.sam
hmh_hu$hm_HNRNPU_flag_hg38_ck_count<-hmh_hu_ck$hm_GFP_flag_hg38_Aligned.out.sam


# calculate rpm by counts*1000000/totalfastqreads and rpkm
hmh_hk$hm_HNRNPK_flag_hg38_merg_rpm<-hmh_hk$hm_HNRNPK_flag_hg38_merg_count*1000000/62978874
hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm<-hmh_hk$hm_HNRNPK_flag_hg38_r1_count*1000000/30746853
hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm<-hmh_hk$hm_HNRNPK_flag_hg38_r2_count*1000000/32232021
hmh_hk$hm_HNRNPK_flag_hg38_ck_rpm<-hmh_hk$hm_HNRNPK_flag_hg38_ck_count*1000000/60319397
hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm_lessCK<-hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm-hmh_hk$hm_HNRNPK_flag_hg38_ck_rpm
hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm_lessCK<-hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm-hmh_hk$hm_HNRNPK_flag_hg38_ck_rpm
hmh_hk$rpkm_r1<-hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm*1000/hmh_hk$Length
hmh_hk$rpkm_r2<-hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm*1000/hmh_hk$Length

hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm_lessCK[which(hmh_hk$hm_HNRNPK_flag_hg38_r1_rpm_lessCK<0)]<-0
hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm_lessCK[which(hmh_hk$hm_HNRNPK_flag_hg38_r2_rpm_lessCK<0)]<-0


hmh_hu$hm_HNRNPU_flag_hg38_merg_rpm<-hmh_hu$hm_HNRNPU_flag_hg38_merg_count*1000000/69381631
hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm<-hmh_hu$hm_HNRNPU_flag_hg38_r1_count*1000000/32401416
hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm<-hmh_hu$hm_HNRNPU_flag_hg38_r2_count*1000000/36980215
hmh_hu$hm_HNRNPU_flag_hg38_ck_rpm<-hmh_hu$hm_HNRNPU_flag_hg38_ck_count*1000000/60319397
hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm_lessCK<-hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm-hmh_hu$hm_HNRNPU_flag_hg38_ck_rpm
hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm_lessCK<-hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm-hmh_hu$hm_HNRNPU_flag_hg38_ck_rpm
hmh_hu$rpkm_r1<-hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm*1000/hmh_hu$Length
hmh_hu$rpkm_r2<-hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm*1000/hmh_hu$Length

hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm_lessCK[which(hmh_hu$hm_HNRNPU_flag_hg38_r1_rpm_lessCK<0)]<-0
hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm_lessCK[which(hmh_hu$hm_HNRNPU_flag_hg38_r2_rpm_lessCK<0)]<-0


write.csv(hmh_hk,'mm10_HNRNPK_hg38peak_10xk_rpkm.csv',row.names = F)
write.csv(hmh_hu,'mm10_HNRNPU_hg38peak_10xk_rpkm.csv',row.names = F)



#################### human
# XIST chrX:73820656-73852723
# KOT1: chr11:2608328-2699994

### HK
(xidx<-which(hmh_hk$Chr=='chrX' & hmh_hk$Strand=='-' & hmh_hk$Start>=73820656 & hmh_hk$End<73852723))
(kidx<-which(hmh_hk$Chr=='chr11' & hmh_hk$Strand=='-' & hmh_hk$Start>=2608328 & hmh_hk$End<2699994))
hmh_hk_10k<-hmh_hk[c(1:10000),]
hmh_hk_x<-hmh_hk[xidx,]
hmh_hk_k<-hmh_hk[kidx,]

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top10k')]<-sum(hmh_hk_10k$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top5k')]<-sum(hmh_hk_10k$rpkm_r1[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top2.5k')]<-sum(hmh_hk_10k$rpkm_r1[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top1k')]<-sum(hmh_hk_10k$rpkm_r1[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='XIST')]<-sum(hmh_hk_x$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='KCNQ1OT1')]<-sum(hmh_hk_k$rpkm_r1)




fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top10k')]<-sum(hmh_hk_10k$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top5k')]<-sum(hmh_hk_10k$rpkm_r2[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top2.5k')]<-sum(hmh_hk_10k$rpkm_r2[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top1k')]<-sum(hmh_hk_10k$rpkm_r2[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='XIST')]<-sum(hmh_hk_x$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPK' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='KCNQ1OT1')]<-sum(hmh_hk_k$rpkm_r2)



### HU
# XIST chrX:73820656-73852723
# KOT1: chr11:2608328-2699994
(xidx<-which(hmh_hu$Chr=='chrX' & hmh_hu$Strand=='-' & hmh_hu$Start>=73820656 & hmh_hu$End<73852723))
(kidx<-which(hmh_hu$Chr=='chr11' & hmh_hu$Strand=='-' & hmh_hu$Start>=2608328 & hmh_hu$End<2699994))
hmh_hu_10k<-hmh_hu[c(1:10000),]
hmh_hu_x<-hmh_hu[xidx,]
hmh_hu_k<-hmh_hu[kidx,]

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top10k')]<-sum(hmh_hu_10k$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top5k')]<-sum(hmh_hu_10k$rpkm_r1[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top2.5k')]<-sum(hmh_hu_10k$rpkm_r1[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='top1k')]<-sum(hmh_hu_10k$rpkm_r1[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='XIST')]<-sum(hmh_hu_x$rpkm_r1)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r1' &
                       fstats$group=='KCNQ1OT1')]<-sum(hmh_hu_k$rpkm_r1)




fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top10k')]<-sum(hmh_hu_10k$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top5k')]<-sum(hmh_hu_10k$rpkm_r2[1:5000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top2.5k')]<-sum(hmh_hu_10k$rpkm_r2[1:2500])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='top1k')]<-sum(hmh_hu_10k$rpkm_r2[1:1000])

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='XIST')]<-sum(hmh_hu_x$rpkm_r2)

fstats$rpkmsum[which(fstats$condition=='HNRNPU' & 
                       fstats$replicate=='human_r2' &
                       fstats$group=='KCNQ1OT1')]<-sum(hmh_hu_k$rpkm_r2)


options(scipen=999)

write.csv(fstats,'HKHU_rpkm_stats.csv',row.names = F)

