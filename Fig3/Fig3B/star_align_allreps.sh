#!/bin/bash

#SBATCH -J peak_%j
#SBATCH -p general
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=12:00:00

## load modules
module load star/2.5.4b
module load samtools


## do following work in directory named $1_classify
mkdir $1_peaks
cd $1_peaks

### 1A. STAR alignment

## STAR alignment of all rips
# move rip files here
mkdir $1
cp ${2} ${1}
gunzip $1/*.fastq.gz

# use STAR aligner to align reads to mm10 mouse genome.
rip_filenames=$(ls -m $1/*.fastq | sed -r 's/\s+//g' | tr -d '\n')
STAR --genomeDir /STAR_genomes/mm10 --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames

# Filter and output to sam file
samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam
