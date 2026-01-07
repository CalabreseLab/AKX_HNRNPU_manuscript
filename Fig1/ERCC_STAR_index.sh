# make star index for ERCC only
module load star/2.7.7a

mkdir ./ERCC_STAR_index

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./ERCC_STAR_index --genomeFastaFiles ERCC92.fa --sjdbGTFfile ERCC92.gtf --genomeSAindexNbases 7
