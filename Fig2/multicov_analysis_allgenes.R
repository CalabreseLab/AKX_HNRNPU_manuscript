# organize multicov results and calculate pearson correlation between proteins

setwd("/work/users/s/h/shuang9/rip/multicov")

library(tidyverse)
library(ggplot2)

header<-read.table('multicov_files_list_allgenes_08092024.txt',sep=' ')
# 68 files in total

header$V1<-gsub('_sorted.bam','',header$V1,fixed=T)

# normalize the reads with total counts in the sequencing fastq files
fastqcount<-read.csv('fastq_read_counts.csv',header=T)
fastqcount$FileName<-gsub('.fastq','',fastqcount$FileName)

# organize fastqcount in the same order of cols in _mc
fastqorder<-header$V1
fastqcount.ordered<-vector(mode='numeric',length=length(fastqorder))


for (n in 1:length(fastqorder)) {
  fname<-fastqorder[n]
  fc<-fastqcount$TotalReads[which(fastqcount$FileName==fname)]
  fastqcount.ordered[n]<-fc
}


# combine replicates of each sample and get mean
sampledata<-read.csv('fastq_samples.csv',header=T)
sampledata$FileName<-gsub('.fastq','',sampledata$FileName)
# only keep the samples are in this analysis
sampledata<-sampledata[which(sampledata$FileName %in% fastqorder),]
uniq_sample<-unique(sampledata$Sample)
# 28 protein in total


lowtriangle<-function(df) {
  
  #df[df<0]<-0
  
  df[upper.tri(df,diag=T)]<-NA
  
  df<-as.data.frame(df)
  
  # add rownames as one column
  df$RBP_1<-rownames(df)
  
  lt<-pivot_longer(df,cols=-RBP_1, names_to='RBP_2', values_to = 'cor_r', values_drop_na = T)
  
  lt<-as.data.frame(lt)
  
  return(lt)
}


mfiles<-dir('/work/users/s/h/shuang9/rip/multicov_out/',pattern='\\.out$')
# 19295 files

# add counter to check if all multicov_out file has 76 columns (6 properties + 70 seq files)
mfile.err<-0

for (n in 1:length(mfiles)) {
  if (n %% 100==0){print(n)}
  #print(n)
  
  fn<-mfiles[n]
  
  f.temp<-read.table(paste0('/work/users/s/h/shuang9/rip/multicov_out/',fn),sep='\t')
  
  # drop V7 and V8 which corresponds to cbx7
  f.temp$V7<-NULL
  f.temp$V8<-NULL
  
  if (ncol(f.temp)!=74) {mfile.err<-mfile.err+1}
  
  colnames(f.temp)<-c('chrom','chromstart','chromend','name','score','strand',header$V1)
  
  # get rpm count
  f_norm<- sweep(f.temp[,c(7:74)], MARGIN=2, STATS=fastqcount.ordered, FUN="/")
  f_norm<- f_norm*1000000
  f_norm<-cbind(f.temp[,c(1:6)],f_norm)
  
  # initialize dataframe to store the mean of each sample
  f_comb<-f_norm[,c(1:6)]
  
  for (sam in uniq_sample) {
    
    fn.temp<-sampledata$FileName[which(sampledata$Sample==sam)]
    
    # only process if there are multiple replicates for a sample
    if (length(fn.temp)>1) {
      
      sam.temp<-f_norm[,which(colnames(f_norm) %in% fn.temp)]
      
      f.mean<-apply(sam.temp,1,mean)
      
      f_comb[[sam]]<-f.mean
      
    } else {
      # just add this sample to the df as is
      f_comb[[sam]]<-f_norm[,which(colnames(f_norm)==fn.temp)]
      
    }
    
  }
  
  fn.comb<-gsub('_multicov.out','_25bp_RBP_comb.csv',fn)
  
  write.csv(f_comb,paste0('/work/users/s/h/shuang9/rip/multicov_RBPcomb/',fn.comb))
  
  # subtract igg from all samples
  # if negative after subtraction, set to 0
  f_comb[,7:34] <- apply(f_comb[,7:34], 2, function(x) x - f_comb$IGG)
  f_comb[,7:34] <- as.data.frame(lapply(f_comb[,7:34], function(x) ifelse(x < 0, 0, x)))
  
  # drop the igg col
  f_comb$IGG<-NULL
  

  # select columns from 7 to 34
  f_data <- f_comb[,7:33]
  
  # compute the pairwise correlations using all complete pairs of observations
  # for each pair of columns, if there are NA entries the result will be NA
  # there is no NA in all datasets
  # NA is produced for RBP that has all 0 counts for a gene for all its cor with others
  f_cor <- cor(f_data, method = "pearson", use = "pairwise.complete.obs")
  
  # apply(f_data, 2, sd)
  
  fn.cor<-gsub('_multicov.out','_25bp_RBP_cor.csv',fn)
  
  write.csv(f_cor,paste0('/work/users/s/h/shuang9/rip/multicov_cor/',fn.cor))
  
  ############ edges 
  
  # Get indices of the upper triangle (excluding the diagonal)
  findices <- which(upper.tri(f_cor,diag = FALSE), arr.ind = TRUE)
  
  # Build data frame
  f_edges <- data.frame(
    Source = rownames(f_cor)[findices[,1]],
    Target = rownames(f_cor)[findices[,2]],
    Weight = abs(f_cor[findices]),
    Color = ifelse(f_cor[findices] > 0, "Positive", "Negative")
  )
  
  
  # Write the edges to a CSV file
  fn.edge<-gsub('_multicov.out','_edges.csv',fn)
  write.csv(f_edges, file = paste0("/work/users/s/h/shuang9/rip/multicov_Gephi/",fn.edge), row.names = FALSE)
  
  ############### generate weighted adjacency matrix for leiden clustering
  
  f_ld<-f_cor
  
  f_ld[f_ld<0]<-0
  
  f_ld[f_ld==1.0]<-0
  
  f_ld[is.na(f_ld)]<-0
  
  fn.ld<-gsub('_multicov.out','_ld_matrix.csv',fn)
  write.csv(f_ld, paste0("/work/users/s/h/shuang9/rip/multicov_Gephi/",fn.ld),row.names = T)
  
  # perform leiden in python as there are warnings
  # of loading the leiden package and warnings of performing leiden
  
  
  ################# integrate all RBP profiles #########################
  # melt the correlation matrix to long format 
  # only keep the upper triangle and exclude the diagonal
  
  # Get indices of the upper triangle (excluding the diagonal)
  findices <- which(upper.tri(f_cor,diag = FALSE), arr.ind = TRUE)
  
  # Build data frame
  f_edges <- data.frame(
    RBP1 = rownames(f_cor)[findices[,1]],
    RBP2 = rownames(f_cor)[findices[,2]],
    Weight = f_cor[findices]
  )

  
  ltname<-gsub('_multicov.out','',fn)
  
  
  if (exists('comb_allgenes')) {
    if(all(comb_allgenes$RBP1==f_edges$RBP1) & all(comb_allgenes$RBP2==f_edges$RBP2)) {
      comb_allgenes[[ltname]]<-f_edges$Weight
    } else {
      print(ltname)
      print('RBP pairs do not match')
    }
    
  } else {
    colnames(f_edges)[3]<-ltname
    comb_allgenes<-f_edges
  }
  
}

write.csv(comb_allgenes,'/work/users/s/h/shuang9/rip/multicov/comb_allgenes_cor_08092024.csv',row.names = F)



RBPpair<-paste0(comb_allgenes$RBP1,'_',comb_allgenes$RBP2)
combdf<-comb_allgenes[,3:19297]
rownames(combdf)<-RBPpair


# keep columns without any NA values
combdf_clean <- combdf[, sapply(combdf, function(col) !any(is.na(col)))]
# 13831

combdf_clean<-as.matrix(combdf_clean)
# If you want to check the rows
zero_var_rows <- apply(combdf_clean, 1, sd) == 0
sum(zero_var_rows) # 0
# If you want to check the columns
zero_var_cols <- apply(combdf_clean, 2, sd) == 0
sum(zero_var_cols) # 0
# Remove rows with zero variation
# combdf <- combdf[!zero_var_rows, ]

# Remove columns with zero variation
combdf.fl <- combdf_clean
rm (combdf_clean)
# 13831


# run leiden cluster with _ld_matrix.csv files in python
# returns _ld_nodes.csv files
# combine all genes ld_nodes info together 
# function to calculate for a certain input list of proteins
# how many genes have them in the same community
genelist<-colnames(combdf.fl)
write.csv(genelist, "genelist13831.csv",row.names = F)

ldfiles<-dir('/work/users/s/h/shuang9/rip/multicov_ld_nodes/')
ldnames<-gsub('_nodes_ld.csv','',ldfiles)
ldfiles<-ldfiles[ldnames %in% genelist]

# a little bit more efficient way using data.table
library(data.table)

ld_allgenes <- NULL

for (n in 1:length(ldfiles)) {
  
  if (n %% 100 == 0) {print(n)}
  
  ldn <- ldfiles[n]
  
  ld.temp <- fread(paste0('/work/users/s/h/shuang9/rip/multicov_ld_nodes/', ldn))
  
  ldgene <- gsub('_nodes_ld.csv', '', ldn)
  
  ld.temp$Label <- NULL
  
  setnames(ld.temp, "Color", ldgene)
  
  if (!is.null(ld_allgenes)) {
    ld_allgenes <- ld_allgenes[ld.temp, on = "Id"]
  } else {
    ld_allgenes <- ld.temp
  }
}

write.csv(ld_allgenes,'/work/users/s/h/shuang9/rip/multicov/ldcommunity_allgenes_08092024.csv',row.names = F)
# 27 x 13831

######################
# iterate thru all possible RBP doublets and count the intra occurrences
samecom<-function(comvec) {
  if (length(unique(comvec))==1) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}

cluster_calculator<-function(proteinvec,ld_allgenes) {
  
  ld.sub<-ld_allgenes[ld_allgenes$Id %in% proteinvec,]
  
  rownames(ld.sub)<-ld.sub$Id
  
  ld.sub$Id<-NULL
  
  ld.stat<-apply(ld.sub,2,samecom)
  
  ld.genelist<-names(ld.stat)[ld.stat==TRUE]
  
  ld.perct<-sum(ld.stat)/length(ld.stat)
  
  ld.count<-sum(ld.stat)
  
  ld.results<-list(ld.perct,ld.count,ld.genelist)
  
  return(ld.results)
  
}


allpro<-ld_allgenes$Id

alldbl<-combn(allpro,2,simplify = F)

alldbl.names<-unlist(lapply(alldbl, function(x) paste(x,collapse='_')))

cluster_loop<-function(namesvec,ld_allgenes) {
  
  res<-as.data.frame(namesvec)
  colnames(res)<-'proset'
  res$perct<-NA
  res$count<-NA
  
  for (n in 1:length(namesvec)) {
    
    if (n %% 100 == 0) {print(n)}
    
    proteinvec<-unlist(strsplit(namesvec[n],'_',fixed=T))
    
    ld.res<-cluster_calculator(proteinvec,ld_allgenes)
    
    res$perct[n]<-ld.res[[1]]
    
    res$count[n]<-ld.res[[2]]
    
  }
  
  return(res)
  
}


alldbl.proset<-cluster_loop(alldbl.names,ld_allgenes)

write.csv(alldbl.proset,'alldoublepro_ld_percent_08092024.csv',row.names = F)
