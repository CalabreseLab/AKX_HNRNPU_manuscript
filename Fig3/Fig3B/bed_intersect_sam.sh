#!/bin/bash

#SBATCH -J ripRank_%j
#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3:00:00


module load samtools

cd $1_peaks

cp /rip/rip_all_stats/comb2igg2reps/$1_peaks_2igg_2reps.bed .

cp /rip/wiggle_2igg2reps/igg_peaks/igg_sorted_neg.bam .
cp /rip/wiggle_2igg2reps/igg_peaks/igg_sorted_pos.bam .
cp /rip/wiggle_2igg2reps/igg_peaks/igg_sorted_neg.bam.bai .
cp /rip/wiggle_2igg2reps/igg_peaks/igg_sorted_pos.bam.bai .

### 1. find aligned reads that overlap with enriched peaks
    # Note: peak file is $1_peaks_2igg.bed

# split peak strands for strand-specific filtering
cat $1_peaks_2igg_2reps.bed | awk '{if ($6=="+") print $0}' > $1_peaks_2igg_pos.bed
cat $1_peaks_2igg_2reps.bed | awk '{if ($6=="-") print $0}' > $1_peaks_2igg_neg.bed

# split sam by strand
samtools view -h -F 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_neg.out.sam
samtools view -h -f 0x10 $1_Aligned_filteredsq30.out.sam > $1_Aligned_pos.out.sam

## convert SAM to BAM
    # -S: take in SAM (by default it expects BAM)
    # -b: output is BAM (by default it produces BAM) 
samtools view -S -b $1_Aligned_neg.out.sam > $1_neg.bam
samtools view -S -b $1_Aligned_pos.out.sam > $1_pos.bam

## sort and index BAM
samtools sort $1_neg.bam -o $1_sorted_neg.bam
samtools sort $1_pos.bam -o $1_sorted_pos.bam
samtools index $1_sorted_neg.bam
samtools index $1_sorted_pos.bam

## filtering: extract SAM entries overlapping peaks.bed from both rip and igg
    # -L: only output alignments overlapping the BED FILE.
    # -h: keep header

samtools view -L "${1}_peaks_2igg_neg.bed" -h -b -o "${1}_peaks_neg.bam" "${1}_sorted_neg.bam"
samtools view -L "${1}_peaks_2igg_pos.bed" -h -b -o "${1}_peaks_pos.bam" "${1}_sorted_pos.bam"

samtools view -L "${1}_peaks_2igg_neg.bed" -h -b -o "igg_${1}_peaks_neg.bam" "igg_sorted_neg.bam"
samtools view -L "${1}_peaks_2igg_pos.bed" -h -b -o "igg_${1}_peaks_pos.bam" "igg_sorted_pos.bam"

scp *_peaks_neg.bam ../bams
scp *_peaks_pos.bam ../bams

