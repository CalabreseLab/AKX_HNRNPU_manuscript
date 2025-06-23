#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

module load bedtools

bedtools getfasta -s -fi /work/users/s/h/shuang9/rip/GRCm38.primary_assembly.genome.fa -bed top1kbedfiles/$1_fc_rpm_2igg_2reps_top1k.bed -fo top1kfasta/$1_fc_rpm_2igg_2reps_top1k.fa