#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBACTH --mem=32g


samtools index $1_sorted.bam
