# concatenate hg38 and mm10 genome together and edit the fasta header to distinguish between h and m

# print fasta headers:
grep '^>' hg38_filtered.fa
grep '^>' mm10_genome_filtered.fa

# add _mm10 to the end of fasta header and concatenate the two fasta files into one
{ sed -e '$a\' hg38_filtered.fa; sed -e '/^>/ s/$/_mm10/' ./mouse/mm10_genome_filtered.fa; } > combined_hg38_mm10.fa
# don't change human chr name as then later it cannot map correctly

grep '^>' combined_hg38_mm10.fa


module load star/2.7.11b
# generate genome indexes
sbatch -p general -t 1-0 -n 24 -N 1 --mem=64g -e staridx.err --wrap='STAR --runMode genomeGenerate --runThreadN 24 --genomeDir STAR_hg38_mm10_index --genomeFastaFiles combined_hg38_mm10.fa'



# change the human and keep the mouse
{ sed -e '$a\' ./mouse/mm10_genome_filtered.fa; sed -e '/^>/ s/$/_hg38/' hg38_filtered.fa; } > combined_mm10_hg38.fa


grep '^>' combined_mm10_hg38.fa


module load star/2.7.11b
# generate genome indexes

sbatch -p general -t 1-0 -n 24 -N 1 --mem=64g -e staridx.err --wrap='STAR --runMode genomeGenerate --runThreadN 24 --genomeDir STAR_mm10_hg38_index --genomeFastaFiles combined_mm10_hg38.fa'




# generate hg38 index by itself
module load star/2.7.11b

sbatch -p general -t 1-0 -n 24 -N 1 --mem=64g -e staridx.err --wrap='STAR --runMode genomeGenerate --runThreadN 24 --genomeDir STAR_hg38_index --genomeFastaFiles hg38_filtered.fa'


