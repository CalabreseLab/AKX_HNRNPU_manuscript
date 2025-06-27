# use Poisson Binomial Distribution to calculate p val
# for protein doublet and triplet
load("/work/users/s/h/shuang9/rip/multicov/allgenes_08092024.RData")

ld_allgenes<-read.csv('ldcommunity_allgenes_08092024.csv',header=T)
# ld_allgenes is 27x13831 that contains leiden community assignments
# for 27 proteins across 13831 genes with the first col as the protein ID

# function to get all the doublet and triplet whether in the same community
# for each gene -- 351 x 13831 or 2925 x 13831
# input is the doublet protein set string
same_com<-function(pro,ld.temp){
  
  # Split the strings in pro and then get the matching shf values for each
  # takes a whole list of protein sets, like dbl$dbl_pro
  procom <- lapply(strsplit(pro, "_", fixed = T), function(x) {
    sapply(x, function(y) ld.temp$ldcom[which(ld.temp$Id == y)])
  })
  
  
  procom_same <- sapply(procom, function(x) {
    if (length(unique(x)) == 1) {
      return(1)
    } else {
      return(0)
    }
  })
  
  return(procom_same)
  
}

allpro<-ld_allgenes$Id
alldbl<-combn(allpro,2,simplify = F)

alldbl.names<-unlist(lapply(alldbl, function(x) paste(x,collapse='_')))

ngene<-ncol(ld_allgenes)-1

# the probability of seeing each protein doublet or triplet 
# being in the same community for a specific gene
alldbl.rate<-vector(mode='numeric',length=ngene)

ld_allgenes<-as.data.frame(ld_allgenes)

for (n in 2:ncol(ld_allgenes)) {
  
  if (n %% 100 == 0) {print(n)}
  
  ld.temp<-ld_allgenes[,c(1,n)]
  
  colnames(ld.temp)[2]<-'ldcom'
  
  dbl.temp<-same_com(alldbl.names,ld.temp)
  alldbl.rate[n-1]<-sum(dbl.temp)/length(dbl.temp)
  
}

write.csv(alldbl.rate,'alldbl_rate_08092024.csv',row.names = F)

#install.packages('poibin')
library(poibin)
library(dplyr)


alldbl.rate<-read.csv('alldbl_rate_08092024.csv',header=T)
alldbl.rate<-alldbl.rate$x

alldbl.proset<-read.csv('alldoublepro_ld_percent_08092024.csv',header=T)

alldbl.proset <- alldbl.proset %>%
  mutate(prev_pval = 1 - ppoibin(count - 1, alldbl.rate),
         rare_pval = ppoibin(count, alldbl.rate))

alldbl.proset$prev_adjp<-p.adjust(alldbl.proset$prev_pval,method='BH')

alldbl.proset$rare_adjp<-p.adjust(alldbl.proset$rare_pval,method='BH')


write.csv(alldbl.proset,'alldoublepro_ld_pval_adj_06182025_raw.csv',row.names = F)



