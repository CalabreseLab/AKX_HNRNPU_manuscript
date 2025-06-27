# generate bedfile for chunks and save score as negative log of the pval

# setwd set working directory properly
setwd("/work/users/s/h/shuang9/rip")
library(Biostrings) # needed to be loaded first, before rtracklayer
library(rtracklayer)
library(data.table) # better data handling



###########################################

# function to take in one row of chunk and exons info and convert to bedfile
convert_transcript_to_genome <- function(temp, exons) {
  
  tempid<-strsplit(temp$labID,'_')
  
  temp_gene<-exons[which(exons$transcript_id==tempid[[1]][7] & exons$gene_id==tempid[[1]][6]),]
  
  if (nrow(temp_gene)==0) {
    
    print('no matched info in gtf')
    print(paste0('transcriptID:',tempid[[1]][7]))
    print(paste0('name:',tempid[[1]][1]))
    break
    
  } else {
    
    sigbed<-temp_gene[,c(1,2,3,5)]
    sigbed$name<-temp$labID
    sigbed$score<-0
    sigbed<-sigbed[,c(1,2,3,5,6,4)]
    colnames(sigbed)<-c('chrom','chromStart','chromEnd','name','score','strand')
    
  }
  
  
  return(sigbed)
  
}


##################################################################
# convert to a bedfile for sig chunks
##################################################################

chunkpval_file<-'TSC-exp_ESC-exp_Chrom-assoc_igg_rpm_4_10_2024_filtered_07222024.csv'
gtf_file<-'gencode.vM25.basic.annotation.gtf'


chunks <- fread(chunkpval_file)
setnames(chunks, "target_id", "labID")

# fix the name with () inside
chunks$labID[grep('Gt(ROSA)26Sor',chunks$labID, fixed=T)]<-'GtROSA26Sor_chr6_113067428_113077333_-_ENSMUSG00000086429.9_ENSMUSG00000086429.9.unspliced'


# extract info from gtf
gtf<-import(gtf_file)

# get the transcript id, gene id, strand, start and end coordinates for each exon
exons <- gtf[ mcols(gtf)$type == "exon" ]

# exons<-as.data.frame(exons)
exons <- as.data.table(exons)

print('finish loading files')


# a function to create BED file based on info from tempbed
# create a BED file for each genes with 25 bp window for that gene
# for each row in tempbed generate a bed file with 25 bp window
# slide from start to end
# name col is: window_number of window
# score col is: 0
# strand col is: strand
# start is the start of the window
# end is the end of the window or the end of the gene, which ever is smaller
# if strand is +, slide from chromStart
# if strand is -, slide from chromEnd
# window_number is the order number of the window
create_bedfile <- function(tempbed, window_size) {
  # Ensure tempbed is a data.table
  tempbed <- as.data.table(tempbed)
  
  # Initialize a data.table to store the new rows
  result <- data.table()
  
  # Initialize window number counter
  window_number <- 1
  
  # Iterate over each row in tempbed
  for (i in 1:nrow(tempbed)) {
    # Extract the current row
    row <- tempbed[i]
    
    # Calculate the start and end positions for each window based on the strand
    if (row$strand == "+") {
      start_pos <- row$chromStart
      end_pos <- row$chromEnd
      while (start_pos < end_pos) {
        current_end <- min(start_pos + window_size - 1, end_pos)
        new_row <- data.table(
          chrom = row$chrom,
          chromStart = start_pos,
          chromEnd = current_end,
          name = paste0(row$name, '_', window_number),
          score = row$score,
          strand = row$strand
        )
        result <- rbind(result, new_row)
        start_pos <- start_pos + window_size
        window_number <- window_number + 1
      }
    } else if (row$strand == "-") {
      start_pos <- row$chromEnd
      end_pos <- row$chromStart
      while (start_pos > end_pos) {
        current_start <- max(start_pos - window_size + 1, end_pos)
        new_row <- data.table(
          chrom = row$chrom,
          chromStart = current_start,
          chromEnd = start_pos,
          name = paste0(row$name, '_', window_number),
          score = row$score,
          strand = row$strand
        )
        result <- rbind(result, new_row)
        start_pos <- start_pos - window_size
        window_number <- window_number + 1
      }
    }
  }
  
  return(result)
}




# set window size 
window_size<-25

# prevent scientific notion of numbers
options(scipen = 999)
# generate bedfiles for all genes in chunks
for (n in 1:nrow(chunks)) {
  
  if (n %% 1000 == 0) {print(n)}
  temp <- chunks[n]
  
  # if gene is unspliced we directly convert the coordinates
  if (grepl('unspliced',temp$labID)) {
    
    tempid<-strsplit(temp$labID,'_')
    
    chrom<-tempid[[1]][2]
    name<-temp$labID
    strand<-tempid[[1]][5]
    score<-0
    chromStart<-as.integer(tempid[[1]][3])
    chromEnd<-as.integer(tempid[[1]][4])
    

    tempbed <- data.table(chrom, chromStart, chromEnd, name, score, strand)
    finaltemp<-create_bedfile(tempbed, window_size)
    fwrite(finaltemp, paste0('bedfiles/',temp$labID,'.bed'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
   
    
  } else {
    # if gene is spliced we need function to calculate coords
    
    conv_temp <- convert_transcript_to_genome(temp, exons) 
    
    tempbed <- as.data.table(conv_temp)
    
    finaltemp<-create_bedfile(tempbed, window_size)
    fwrite(finaltemp, paste0('bedfiles/',temp$labID,'.bed'), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    

    
  }
  
}







