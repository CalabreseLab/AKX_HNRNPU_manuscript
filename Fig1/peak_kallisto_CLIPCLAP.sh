#!/bin/bash

#SBATCH -p general
#SBATCH --mem=120G
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=50:00:00


module load kallisto/0.51.1
module load samtools/1.22


mkdir CLIPCLAP_kallisto
cd CLIPCLAP_kallisto


# copy over the true peaks bedfiles
cp ../SAFA_MinusTag_CLAP_mm10_truepeaks.bed .
cp ../SAFA_PlusTag_CLAP_hg38_truepeaks.bed .

### 1. find aligned reads that overlap with enriched peaks

# split peak strands for strand-specific filtering
cat SAFA_MinusTag_CLAP_mm10_truepeaks.bed | awk '{if ($6=="+") print $0}' > SAFA_MinusTag_CLAP_mm10_truepeaks_pos.bed
cat SAFA_MinusTag_CLAP_mm10_truepeaks.bed | awk '{if ($6=="-") print $0}' > SAFA_MinusTag_CLAP_mm10_truepeaks_neg.bed

cat SAFA_PlusTag_CLAP_hg38_truepeaks.bed | awk '{if ($6=="+") print $0}' > SAFA_PlusTag_CLAP_hg38_truepeaks_pos.bed
cat SAFA_PlusTag_CLAP_hg38_truepeaks.bed | awk '{if ($6=="-") print $0}' > SAFA_PlusTag_CLAP_hg38_truepeaks_neg.bed


## sort and index BAM
samtools sort ../strands/SAFA_MinusTag_CLAP_mm10_pairedonly_neg.bam -o SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_MinusTag_CLAP_mm10_pairedonly_pos.bam -o SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_pos.bam
samtools index SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_neg.bam
samtools index SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_pos.bam

samtools sort ../strands/SAFA_PlusTag_CLAP_mm10_pairedonly_neg.bam -o SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_PlusTag_CLAP_mm10_pairedonly_pos.bam -o SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_pos.bam
samtools index SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_neg.bam
samtools index SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_pos.bam



samtools sort ../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_neg.bam -o SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_pos.bam -o SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_pos.bam
samtools index SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_neg.bam
samtools index SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_pos.bam

samtools sort ../strands/SAFA_MinusTag_CLAP_hg38_pairedonly_neg.bam -o SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_MinusTag_CLAP_hg38_pairedonly_pos.bam -o SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_pos.bam
samtools index SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_neg.bam
samtools index SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_pos.bam




samtools sort ../strands/SAFA_MinusTag_CLIP_mm10_pairedonly_neg.bam -o SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_MinusTag_CLIP_mm10_pairedonly_pos.bam -o SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_pos.bam
samtools index SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_neg.bam
samtools index SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_pos.bam

samtools sort ../strands/SAFA_PlusTag_CLIP_mm10_pairedonly_neg.bam -o SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_PlusTag_CLIP_mm10_pairedonly_pos.bam -o SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_pos.bam
samtools index SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_neg.bam
samtools index SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_pos.bam



samtools sort ../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_neg.bam -o SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_pos.bam -o SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_pos.bam
samtools index SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_neg.bam
samtools index SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_pos.bam

samtools sort ../strands/SAFA_MinusTag_CLIP_hg38_pairedonly_neg.bam -o SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_neg.bam
samtools sort ../strands/SAFA_MinusTag_CLIP_hg38_pairedonly_pos.bam -o SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_pos.bam
samtools index SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_neg.bam
samtools index SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_pos.bam



## filtering: extract SAM entries overlapping peaks.bed from both rip and igg
    # -L: only output alignments overlapping the BED FILE.
    # -h: keep header
    # -b: output bam file
    # --fetch-pairs ensure one read overlaps my target region; still output its mate even if the mate is outside the region
    # -f 2 Keep only properly paired reads for the output bam
samtools view SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_neg.bed > SAFA_MinusTag_CLAP_mm10_peaks_neg.out.bam
samtools view SAFA_MinusTag_CLAP_mm10_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_pos.bed > SAFA_MinusTag_CLAP_mm10_peaks_pos.out.bam

samtools view SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_neg.bed > SAFA_MinusTag_CLAP_mm10_CKpeaks_neg.out.bam
samtools view SAFA_PlusTag_CLAP_mm10_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_pos.bed > SAFA_MinusTag_CLAP_mm10_CKpeaks_pos.out.bam



samtools view SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_neg.bed > SAFA_PlusTag_CLAP_hg38_peaks_neg.out.bam
samtools view SAFA_PlusTag_CLAP_hg38_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_pos.bed > SAFA_PlusTag_CLAP_hg38_peaks_pos.out.bam

samtools view SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_neg.bed > SAFA_PlusTag_CLAP_hg38_CKpeaks_neg.out.bam
samtools view SAFA_MinusTag_CLAP_hg38_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_pos.bed > SAFA_PlusTag_CLAP_hg38_CKpeaks_pos.out.bam




samtools view SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_neg.bed > SAFA_MinusTag_CLIP_mm10_peaks_neg.out.bam
samtools view SAFA_MinusTag_CLIP_mm10_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_pos.bed > SAFA_MinusTag_CLIP_mm10_peaks_pos.out.bam

samtools view SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_neg.bed > SAFA_MinusTag_CLIP_mm10_CKpeaks_neg.out.bam
samtools view SAFA_PlusTag_CLIP_mm10_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_MinusTag_CLAP_mm10_truepeaks_pos.bed > SAFA_MinusTag_CLIP_mm10_CKpeaks_pos.out.bam



samtools view SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_neg.bed > SAFA_PlusTag_CLIP_hg38_peaks_neg.out.bam
samtools view SAFA_PlusTag_CLIP_hg38_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_pos.bed > SAFA_PlusTag_CLIP_hg38_peaks_pos.out.bam

samtools view SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_neg.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_neg.bed > SAFA_PlusTag_CLIP_hg38_CKpeaks_neg.out.bam
samtools view SAFA_MinusTag_CLIP_hg38_pairedonly_sorted_pos.bam -h -b --fetch-pairs -f 2 -L SAFA_PlusTag_CLAP_hg38_truepeaks_pos.bed > SAFA_PlusTag_CLIP_hg38_CKpeaks_pos.out.bam



# need to merge the pos and neg bam into one bam before convert to fastq
samtools merge -@ 8 SAFA_MinusTag_CLAP_mm10_peaks_merged.bam SAFA_MinusTag_CLAP_mm10_peaks_pos.out.bam SAFA_MinusTag_CLAP_mm10_peaks_neg.out.bam
samtools merge -@ 8 SAFA_MinusTag_CLAP_mm10_CKpeaks_merged.bam SAFA_MinusTag_CLAP_mm10_CKpeaks_pos.out.bam SAFA_MinusTag_CLAP_mm10_CKpeaks_neg.out.bam
# sort before conversion
samtools sort -@ 8 -n -o SAFA_MinusTag_CLAP_mm10_peaks_merged_sorted.bam SAFA_MinusTag_CLAP_mm10_peaks_merged.bam
samtools sort -@ 8 -n -o SAFA_MinusTag_CLAP_mm10_CKpeaks_merged_sorted.bam SAFA_MinusTag_CLAP_mm10_CKpeaks_merged.bam



samtools merge -@ 8 SAFA_PlusTag_CLAP_hg38_peaks_merged.bam SAFA_PlusTag_CLAP_hg38_peaks_pos.out.bam SAFA_PlusTag_CLAP_hg38_peaks_neg.out.bam
samtools merge -@ 8 SAFA_PlusTag_CLAP_hg38_CKpeaks_merged.bam SAFA_PlusTag_CLAP_hg38_CKpeaks_pos.out.bam SAFA_PlusTag_CLAP_hg38_CKpeaks_neg.out.bam
# sort before conversion
samtools sort -@ 8 -n -o SAFA_PlusTag_CLAP_hg38_peaks_merged_sorted.bam SAFA_PlusTag_CLAP_hg38_peaks_merged.bam
samtools sort -@ 8 -n -o SAFA_PlusTag_CLAP_hg38_CKpeaks_merged_sorted.bam SAFA_PlusTag_CLAP_hg38_CKpeaks_merged.bam





samtools merge -@ 8 SAFA_MinusTag_CLIP_mm10_peaks_merged.bam SAFA_MinusTag_CLIP_mm10_peaks_pos.out.bam SAFA_MinusTag_CLIP_mm10_peaks_neg.out.bam
samtools merge -@ 8 SAFA_MinusTag_CLIP_mm10_CKpeaks_merged.bam SAFA_MinusTag_CLIP_mm10_CKpeaks_pos.out.bam SAFA_MinusTag_CLIP_mm10_CKpeaks_neg.out.bam
# sort before conversion
samtools sort -@ 8 -n -o SAFA_MinusTag_CLIP_mm10_peaks_merged_sorted.bam SAFA_MinusTag_CLIP_mm10_peaks_merged.bam
samtools sort -@ 8 -n -o SAFA_MinusTag_CLIP_mm10_CKpeaks_merged_sorted.bam SAFA_MinusTag_CLIP_mm10_CKpeaks_merged.bam



samtools merge -@ 8 SAFA_PlusTag_CLIP_hg38_peaks_merged.bam SAFA_PlusTag_CLIP_hg38_peaks_pos.out.bam SAFA_PlusTag_CLIP_hg38_peaks_neg.out.bam
samtools merge -@ 8 SAFA_PlusTag_CLIP_hg38_CKpeaks_merged.bam SAFA_PlusTag_CLIP_hg38_CKpeaks_pos.out.bam SAFA_PlusTag_CLIP_hg38_CKpeaks_neg.out.bam
# sort before conversion
samtools sort -@ 8 -n -o SAFA_PlusTag_CLIP_hg38_peaks_merged_sorted.bam SAFA_PlusTag_CLIP_hg38_peaks_merged.bam
samtools sort -@ 8 -n -o SAFA_PlusTag_CLIP_hg38_CKpeaks_merged_sorted.bam SAFA_PlusTag_CLIP_hg38_CKpeaks_merged.bam



# convert BAM to FASTQ
# but here we need to generate two fastq as the samples are paired end
# -0 /dev/null drops reads that aren’t in a proper pair as “singletons”
# -s /dev/null drops orphaned reads
# -n keeps read names consistent
samtools fastq -@ 8 -n -1 SAFA_MinusTag_CLAP_mm10_peaks_R1.fastq -2 SAFA_MinusTag_CLAP_mm10_peaks_R2.fastq -0 /dev/null -s /dev/null SAFA_MinusTag_CLAP_mm10_peaks_merged_sorted.bam
samtools fastq -@ 8 -n -1 SAFA_MinusTag_CLAP_mm10_CKpeaks_R1.fastq -2 SAFA_MinusTag_CLAP_mm10_CKpeaks_R2.fastq -0 /dev/null -s /dev/null SAFA_MinusTag_CLAP_mm10_CKpeaks_merged_sorted.bam


samtools fastq -@ 8 -n -1 SAFA_PlusTag_CLAP_hg38_peaks_R1.fastq -2 SAFA_PlusTag_CLAP_hg38_peaks_R2.fastq -0 /dev/null -s /dev/null SAFA_PlusTag_CLAP_hg38_peaks_merged_sorted.bam
samtools fastq -@ 8 -n -1 SAFA_PlusTag_CLAP_hg38_CKpeaks_R1.fastq -2 SAFA_PlusTag_CLAP_hg38_CKpeaks_R2.fastq -0 /dev/null -s /dev/null SAFA_PlusTag_CLAP_hg38_CKpeaks_merged_sorted.bam



samtools fastq -@ 8 -n -1 SAFA_MinusTag_CLIP_mm10_peaks_R1.fastq -2 SAFA_MinusTag_CLIP_mm10_peaks_R2.fastq -0 /dev/null -s /dev/null SAFA_MinusTag_CLIP_mm10_peaks_merged_sorted.bam
samtools fastq -@ 8 -n -1 SAFA_MinusTag_CLIP_mm10_CKpeaks_R1.fastq -2 SAFA_MinusTag_CLIP_mm10_CKpeaks_R2.fastq -0 /dev/null -s /dev/null SAFA_MinusTag_CLIP_mm10_CKpeaks_merged_sorted.bam


samtools fastq -@ 8 -n -1 SAFA_PlusTag_CLIP_hg38_peaks_R1.fastq -2 SAFA_PlusTag_CLIP_hg38_peaks_R2.fastq -0 /dev/null -s /dev/null SAFA_PlusTag_CLIP_hg38_peaks_merged_sorted.bam
samtools fastq -@ 8 -n -1 SAFA_PlusTag_CLIP_hg38_CKpeaks_R1.fastq -2 SAFA_PlusTag_CLIP_hg38_CKpeaks_R2.fastq -0 /dev/null -s /dev/null SAFA_PlusTag_CLIP_hg38_CKpeaks_merged_sorted.bam


### 2. Kallisto: align both rip and igg peak.fastq to transcriptome

kallisto_idx_m='kallisto_transcriptome_vM25_01132026.index'
kallisto_idx_h='kallisto_transcriptome_v47_01132026.index'

## kallisto alignment
    # -i: location of the transcriptome index
    # -o: kallisto output folder for each dataset
    # --single -l 200 -s 50: single-end data
    # -l and -s are required for SE reads: estimated average and standard deviation of fragment length
    # rf-stranded: reverse stranded
    # -b 100 quantify uncertainty using 100 bootstrap resamples
kallisto quant -i $kallisto_idx_m -o SAFA_MinusTag_CLAP_mm10_peaks_kallisto --rf-stranded SAFA_MinusTag_CLAP_mm10_peaks_R1.fastq SAFA_MinusTag_CLAP_mm10_peaks_R2.fastq
kallisto quant -i $kallisto_idx_m -o SAFA_MinusTag_CLAP_mm10_CKpeaks_kallisto --rf-stranded SAFA_MinusTag_CLAP_mm10_CKpeaks_R1.fastq SAFA_MinusTag_CLAP_mm10_CKpeaks_R2.fastq

kallisto quant -i $kallisto_idx_h -o SAFA_PlusTag_CLAP_hg38_peaks_kallisto --rf-stranded SAFA_PlusTag_CLAP_hg38_peaks_R1.fastq SAFA_PlusTag_CLAP_hg38_peaks_R2.fastq
kallisto quant -i $kallisto_idx_h -o SAFA_PlusTag_CLAP_hg38_CKpeaks_kallisto --rf-stranded SAFA_PlusTag_CLAP_hg38_CKpeaks_R1.fastq SAFA_PlusTag_CLAP_hg38_CKpeaks_R2.fastq



kallisto quant -i $kallisto_idx_m -o SAFA_MinusTag_CLIP_mm10_peaks_kallisto --rf-stranded SAFA_MinusTag_CLIP_mm10_peaks_R1.fastq SAFA_MinusTag_CLIP_mm10_peaks_R2.fastq
kallisto quant -i $kallisto_idx_m -o SAFA_MinusTag_CLIP_mm10_CKpeaks_kallisto --rf-stranded SAFA_MinusTag_CLIP_mm10_CKpeaks_R1.fastq SAFA_MinusTag_CLIP_mm10_CKpeaks_R2.fastq

kallisto quant -i $kallisto_idx_g -o SAFA_PlusTag_CLIP_hg38_peaks_kallisto --rf-stranded SAFA_PlusTag_CLIP_hg38_peaks_R1.fastq SAFA_PlusTag_CLIP_hg38_peaks_R2.fastq
kallisto quant -i $kallisto_idx_g -o SAFA_PlusTag_CLIP_hg38_CKpeaks_kallisto --rf-stranded SAFA_PlusTag_CLIP_hg38_CKpeaks_R1.fastq SAFA_PlusTag_CLIP_hg38_CKpeaks_R2.fastq



