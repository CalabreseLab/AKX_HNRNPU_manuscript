# get the top 1000 peaks of each RBP after -igg
# generate bedfiles for the peaks

# get all files under comb2igg2reps
filelist<-list.files('comb2igg2reps/')

# loop thru each RBP
for (f in filelist) {
  
  pname<-strsplit(f,'_',fixed=T)[[1]][1]

  p_rpm<-paste0(pname,'_rpm')
  
  df<-read.table(paste0('comb2igg2reps/',f),sep='\t',header=T)
  
  print(which(colnames(df)==p_rpm))
  
  df$overigg<-df[[p_rpm]]-df$igg_rpm
  
  df_sorted <- df[order(-df$overigg), ]
  
  df_1k<-df_sorted[c(1:1000),]
  
  df_1k$score<-0
  
  # generate bedfile 
  
  bed<-cbind(df_1k$Chr,df_1k$Start,df_1k$End,df_1k$GeneID,df_1k$score,df_1k$Strand)
  
  bedname<-gsub('.txt','_top1k.bed',f)
  
  write.table(bed,paste0('top1kbedfiles/',bedname),sep='\t',row.names = F,col.names = F, quote=F)
  
}
