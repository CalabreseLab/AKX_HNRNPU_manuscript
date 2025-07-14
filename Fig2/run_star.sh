#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=64g


STAR --genomeDir /rip/MM10_genome_index --runThreadN 8 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesIn /rip/rip_data/$1.fastq --outFileNamePrefix /rip/star_aligned/$1_
