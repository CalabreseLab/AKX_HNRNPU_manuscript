#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

samtools sort $1.bam -o $1_sorted.bam
