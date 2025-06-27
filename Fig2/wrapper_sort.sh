#!/bin/bash
#SBATCH --job-name=bam_sort
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

module load samtools

for bam_file in *.bam
do
  base=${bam_file%.bam}
  sbatch --export=ALL,base=${base} sort_bam.sh ${base}
done
