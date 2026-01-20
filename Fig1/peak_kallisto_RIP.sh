#!/bin/bash

#SBATCH -p general
#SBATCH --mem=120G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=20:00:00


module load kallisto/0.51.1
module load samtools/1.22


mkdir RIP_kallisto
cd RIP_kallisto

# copy over the true peaks bedfile
cp ../hm_HNRNPK_flag_truepeaks.bed .
cp ../hm_HNRNPU_flag_truepeaks.bed .

### 1. find aligned reads that overlap with true peaks

# split peak strands for strand-specific filtering
cat hm_HNRNPK_flag_truepeaks.bed | awk '{if ($6=="+") print $0}' > hm_HNRNPK_flag_truepeaks_pos.bed
cat hm_HNRNPK_flag_truepeaks.bed | awk '{if ($6=="-") print $0}' > hm_HNRNPK_flag_truepeaks_neg.bed

cat hm_HNRNPU_flag_truepeaks.bed | awk '{if ($6=="+") print $0}' > hm_HNRNPU_flag_truepeaks_pos.bed
cat hm_HNRNPU_flag_truepeaks.bed | awk '{if ($6=="-") print $0}' > hm_HNRNPU_flag_truepeaks_neg.bed


## sort and index BAM
samtools sort ../strands/hm_GFP_flag_Aligned.out_neg.bam -o hm_GFP_flag_sorted_neg.bam
samtools sort ../strands/hm_GFP_flag_Aligned.out_pos.bam -o hm_GFP_flag_sorted_pos.bam
samtools index hm_GFP_flag_sorted_neg.bam
samtools index hm_GFP_flag_sorted_pos.bam


samtools sort ../strands/hm_HNRNPK_flag_Aligned.out_neg.bam -o hm_HNRNPK_flag_sorted_neg.bam
samtools sort ../strands/hm_HNRNPK_flag_Aligned.out_pos.bam -o hm_HNRNPK_flag_sorted_pos.bam
samtools index hm_HNRNPK_flag_sorted_neg.bam
samtools index hm_HNRNPK_flag_sorted_pos.bam


samtools sort ../strands/hm_HNRNPU_flag_Aligned.out_neg.bam -o hm_HNRNPU_flag_sorted_neg.bam
samtools sort ../strands/hm_HNRNPU_flag_Aligned.out_pos.bam -o hm_HNRNPU_flag_sorted_pos.bam
samtools index hm_HNRNPU_flag_sorted_neg.bam
samtools index hm_HNRNPU_flag_sorted_pos.bam


## filtering: extract SAM entries overlapping peaks.bed from both rip and control
    # -L: only output alignments overlapping the BED FILE.
    # -h: keep header
samtools view hm_HNRNPK_flag_sorted_neg.bam -h -b -L hm_HNRNPK_flag_truepeaks_neg.bed > hm_HNRNPK_flag_peaks_neg.bam
samtools view hm_HNRNPK_flag_sorted_pos.bam -h -b -L hm_HNRNPK_flag_truepeaks_pos.bed > hm_HNRNPK_flag_peaks_pos.bam

samtools view hm_GFP_flag_sorted_neg.bam -h -b -L hm_HNRNPK_flag_truepeaks_neg.bed > hm_GFP_flag_HKpeaks_neg.bam
samtools view hm_GFP_flag_sorted_pos.bam -h -b -L hm_HNRNPK_flag_truepeaks_pos.bed > hm_GFP_flag_HKpeaks_pos.bam



samtools view hm_HNRNPU_flag_sorted_neg.bam -h -b -L hm_HNRNPU_flag_truepeaks_neg.bed > hm_HNRNPU_flag_peaks_neg.bam
samtools view hm_HNRNPU_flag_sorted_pos.bam -h -b -L hm_HNRNPU_flag_truepeaks_pos.bed > hm_HNRNPU_flag_peaks_pos.bam

samtools view hm_GFP_flag_sorted_neg.bam -h -b -L hm_HNRNPU_flag_truepeaks_neg.bed > hm_GFP_flag_HUpeaks_neg.bam
samtools view hm_GFP_flag_sorted_pos.bam -h -b -L hm_HNRNPU_flag_truepeaks_pos.bed > hm_GFP_flag_HUpeaks_pos.bam



# need to merge the pos and neg bam into one bam before convert to fastq
samtools merge -@ 8 hm_HNRNPK_flag_peaks_merged.bam hm_HNRNPK_flag_peaks_pos.bam hm_HNRNPK_flag_peaks_neg.bam
samtools merge -@ 8 hm_GFP_flag_HKpeaks_merged.bam hm_GFP_flag_HKpeaks_pos.bam hm_GFP_flag_HKpeaks_neg.bam
# sort before conversion
samtools sort -@ 8 -n -o hm_HNRNPK_flag_peaks_merged_sorted.bam hm_HNRNPK_flag_peaks_merged.bam
samtools sort -@ 8 -n -o hm_GFP_flag_HKpeaks_merged_sorted.bam hm_GFP_flag_HKpeaks_merged.bam



# need to merge the pos and neg bam into one bam before convert to fastq
samtools merge -@ 8 hm_HNRNPU_flag_peaks_merged.bam hm_HNRNPU_flag_peaks_pos.bam hm_HNRNPU_flag_peaks_neg.bam
samtools merge -@ 8 hm_GFP_flag_HUpeaks_merged.bam hm_GFP_flag_HUpeaks_pos.bam hm_GFP_flag_HUpeaks_neg.bam
# sort before conversion
samtools sort -@ 8 -n -o hm_HNRNPU_flag_peaks_merged_sorted.bam hm_HNRNPU_flag_peaks_merged.bam
samtools sort -@ 8 -n -o hm_GFP_flag_HUpeaks_merged_sorted.bam hm_GFP_flag_HUpeaks_merged.bam


# convert BAM to FASTQ
samtools fastq hm_HNRNPK_flag_peaks_merged_sorted.bam > hm_HNRNPK_flag_peaks.fastq

samtools fastq hm_GFP_flag_HKpeaks_merged_sorted.bam > hm_GFP_flag_HKpeaks.fastq

samtools fastq hm_HNRNPU_flag_peaks_merged_sorted.bam > hm_HNRNPU_flag_peaks.fastq

samtools fastq hm_GFP_flag_HUpeaks_merged_sorted.bam > hm_GFP_flag_HUpeaks.fastq


### 2. Kallisto: align both rip and control peak.fastq to transcriptome

kallisto_idx='kallisto_transcriptome_vM25_01132026.index'


## kallisto alignment
    # -i: location of the transcriptome index
    # -o: kallisto output folder for each dataset
    # --single -l 200 -s 50: single-end data
    # -l and -s are required for SE reads: estimated average and standard deviation of fragment length
    # rf-stranded: reverse stranded
kallisto quant -i $kallisto_idx -o hm_HNRNPK_kallisto --single -l 200 -s 50 --rf-stranded hm_HNRNPK_flag_peaks.fastq
kallisto quant -i $kallisto_idx -o hm_HNRNPK_CK_kallisto --single -l 200 -s 50 --rf-stranded hm_GFP_flag_HKpeaks.fastq

kallisto quant -i $kallisto_idx -o hm_HNRNPU_kallisto --single -l 200 -s 50 --rf-stranded hm_HNRNPU_flag_peaks.fastq
kallisto quant -i $kallisto_idx -o hm_HNRNPU_CK_kallisto --single -l 200 -s 50 --rf-stranded hm_GFP_flag_HUpeaks.fastq

