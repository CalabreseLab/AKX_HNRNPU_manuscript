#!/bin/bash

#SBATCH -J igg
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3:00:00

module load samtools

cd igg_peaks


### generate igg bam split by strand

# split sam by strand
samtools view -h -F 0x10 igg_Aligned_filteredsq30.out.sam > igg_Aligned_neg.out.sam
samtools view -h -f 0x10 igg_Aligned_filteredsq30.out.sam > igg_Aligned_pos.out.sam

## convert SAM to BAM
    # -S: take in SAM (by default it expects BAM)
    # -b: output is BAM (by default it produces BAM) 
samtools view -S -b igg_Aligned_neg.out.sam > igg_neg.bam
samtools view -S -b igg_Aligned_pos.out.sam > igg_pos.bam

## sort and index BAM
samtools sort igg_neg.bam -o igg_sorted_neg.bam
samtools sort igg_pos.bam -o igg_sorted_pos.bam
samtools index igg_sorted_neg.bam
samtools index igg_sorted_pos.bam

