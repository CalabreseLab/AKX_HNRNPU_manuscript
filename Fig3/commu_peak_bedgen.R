
library(data.table)
library(dplyr)

######## Xist
# Xist chrX:103,460,366-103,483,254 -

files<-dir(pattern='^Xist_.*\\.wig$')

for (f in files) {
  
  xtemp<-fread(f,sep='\t',skip=1)
  xtemp$name<-''
  xtemp$strand<-'-'
  xtemp<-xtemp[,c(1:3,5,4,6)]
  xtemp<-xtemp[V4>0.125]
  
  xname<-gsub('wig','bed',f)
  
  fwrite(xtemp,xname,sep='\t',quote=F,row.names = F,col.names = F)

}

######## Airn
# Airn: chr17:12,741,398-12,830,151 +
files<-dir(pattern='^Airn_.*\\.wig$')

for (f in files) {
  
  xtemp<-fread(f,sep='\t',skip=1)
  xtemp$name<-''
  xtemp$strand<-'+'
  xtemp<-xtemp[,c(1:3,5,4,6)]
  xtemp<-xtemp[V4>0.125]

  xname<-gsub('wig','bed',f)
  
  fwrite(xtemp,xname,sep='\t',quote=F,row.names = F,col.names = F)
  
}

######## Kot1
# Kot1: chr7:143,203,458-143,296,549 -
# exclude c6 as it is empty

files<-dir(pattern='^Kot1_.*\\.wig$')
files<-files[-7] # remove c6
# order it properly
files<-files[c(1,3:9,2)]

for (f in files) {
  
  xtemp<-fread(f,sep='\t',skip=1)
  xtemp$name<-''
  xtemp$strand<-'-'
  xtemp<-xtemp[,c(1:3,5,4,6)]
  xtemp<-xtemp[V4>0.125]
  
  if (nrow(xtemp)>0) {
    xname<-gsub('wig','bed',f)
    
    fwrite(xtemp,xname,sep='\t',quote=F,row.names = F,col.names = F)
    
  }
  
  
}

