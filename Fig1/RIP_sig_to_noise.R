# for hm_ samples
# get signal to noise bar plot for X and K in Flag samples 
# for HK and HU calculate rpmlessCK and CK and plot barplot
# sum all rpms under true peaks that falls in X and K

mm10_hk<-read.csv('mm10_HNRNPK_all_rpkm.csv',header=T)
mm10_hu<-read.csv('mm10_HNRNPU_all_rpkm.csv',header=T)


# human
# XIST chrX:73820656-73852723
# KOT1: chr11:2608328-2699994

get_hx<-function(df) {
  xidx<-which(df$Chr=='chrX' & df$Strand=='-' & df$Start>=73820656 & df$End<73852723)
  return(df[xidx,])
}

get_hot1<-function(df){
  kidx<-which(df$Chr=='chr11' & df$Strand=='-' & df$Start>=2608328 & df$End<2699994)
  return(df[kidx,])
}

# mouse
# Xist: chrX:103460366-103483254
# Kot1: chr7:143203458-143296549

get_mx<-function(df){
  xidx<-which(df$Chr=='chrX' & df$Strand=='-' & df$Start>=103460366 & df$End<103483254)
  return(df[xidx,])
}

get_mot1<-function(df){
  kidx<-which(df$Chr=='chr7' & df$Strand=='-' & df$Start>=143203458 & df$End<143296549)
  return(df[kidx,])
}


# get all stats for xist
xstats<-matrix(nrow=6,ncol=4)
xstats<-as.data.frame(xstats)
colnames(xstats)<-c('condition','replicate','groups','rpmsum')
xstats$condition<-rep(c('mm10_HNRNPK','mm10_HNRNPU'),each=3)
xstats$replicate<-rep(c('r1','r2','r1'),times=2)
xstats$groups<-rep(c('RBPlessCK','RBPlessCK','CK'),times=2)


mm_hnrnpk<-get_mx(mm10_hk)
mm_hnrnpu<-get_mx(mm10_hu)

  
xstats$rpmsum[1]<-sum(mm_hnrnpk$hm_HNRNPK_flag_r1_rpm_lessCK)
xstats$rpmsum[2]<-sum(mm_hnrnpk$hm_HNRNPK_flag_r2_rpm_lessCK)
xstats$rpmsum[3]<-sum(mm_hnrnpk$hm_HNRNPK_flag_ck_rpm)

xstats$rpmsum[4]<-sum(mm_hnrnpu$hm_HNRNPU_flag_r1_rpm_lessCK)
xstats$rpmsum[5]<-sum(mm_hnrnpu$hm_HNRNPU_flag_r2_rpm_lessCK)
xstats$rpmsum[6]<-sum(mm_hnrnpu$hm_HNRNPU_flag_ck_rpm)

write.csv(xstats,'Xist_RIP_signaltonoise_points.csv',row.names = F)


# get all stats for kot1
xstats<-matrix(nrow=6,ncol=4)
xstats<-as.data.frame(xstats)
colnames(xstats)<-c('condition','replicate','groups','rpmsum')
xstats$condition<-rep(c('mm10_HNRNPK','mm10_HNRNPU'),each=3)
xstats$replicate<-rep(c('r1','r2','r1'),times=2)
xstats$groups<-rep(c('RBPlessCK','RBPlessCK','CK'),times=2)


mm_hnrnpk<-get_mot1(mm10_hk)
mm_hnrnpu<-get_mot1(mm10_hu)

xstats$rpmsum[1]<-sum(mm_hnrnpk$hm_HNRNPK_flag_r1_rpm_lessCK)
xstats$rpmsum[2]<-sum(mm_hnrnpk$hm_HNRNPK_flag_r2_rpm_lessCK)
xstats$rpmsum[3]<-sum(mm_hnrnpk$hm_HNRNPK_flag_ck_rpm)

xstats$rpmsum[4]<-sum(mm_hnrnpu$hm_HNRNPU_flag_r1_rpm_lessCK)
xstats$rpmsum[5]<-sum(mm_hnrnpu$hm_HNRNPU_flag_r2_rpm_lessCK)
xstats$rpmsum[6]<-sum(mm_hnrnpu$hm_HNRNPU_flag_ck_rpm)



write.csv(xstats,'Kot1_RIP_signaltonoise_points.csv',row.names = F)

# merge with CLIP CLAP results and plot over there
