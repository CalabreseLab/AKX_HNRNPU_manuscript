#!/bin/bash
#SBATCH --job-name=bam_index
#SBATCH --output=bam_index_%j.out
#SBATCH --error=bam_index_%j.err
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBACTH --mem=32g

module load samtools

for bam_file in *_sorted.bam
do
  base=${bam_file%_sorted.bam}
  sbatch --export=ALL,base=${base} index_bam.sh ${base}
done
