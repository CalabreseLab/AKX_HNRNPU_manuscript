#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

module load meme

streme --oc top1kmotif/streme_$1 -n top1kcontrol/control_$1_fc_rpm_2igg_2reps_top1k.fa --maxw 8 --minw 4 --nmotifs 10 --rna --p top1kfasta/$1_fc_rpm_2igg_2reps_top1k.fa