#!/bin/bash

## load modules
module load star/2.7.11b
module load samtools/1.22
module load macs/2.2.9.1
module load subread/2.1.1
module load bedtools/2.31.1


#### Part 1: from rip, STAR genome alignment and MACS peak-calling 

### 1A. STAR alignment

## STAR alignment of all rips

# samples that has both human and mouse cells
# map to mouse

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_r1_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r1_rip_5_S5_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_r2_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r2_rip_6_S6_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_r1_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r1_rip_5_S11_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_r2_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r2_rip_6_S12_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_r1_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r1_rip_7_S7_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_r2_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r2_rip_8_S8_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_r1_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r1_rip_7_S13_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_r2_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r2_rip_8_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_r1_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r1_rip_9_S9_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_r2_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r2_rip_10_S10_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_r1_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r1_rip_9_S15_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_r2_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r2_rip_10_S16_R1_001.fastq"


# samples that has both human and mouse cells
# map to human

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_r1_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r1_rip_5_S5_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_r2_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r2_rip_6_S6_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_r1_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r1_rip_5_S11_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_r2_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r2_rip_6_S12_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_r1_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r1_rip_7_S7_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_r2_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r2_rip_8_S8_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_r1_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r1_rip_7_S13_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_r2_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r2_rip_8_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_r1_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r1_rip_9_S9_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_r2_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r2_rip_10_S10_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_r1_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r1_rip_9_S15_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_r2_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r2_rip_10_S16_R1_001.fastq"


# samples that has only human cells
# mapped to human only index

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPU_r1_ --readFilesIn G10_293T_hnrnpu_r1_rip_1_S1_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPU_r2_ --readFilesIn G10_293T_hnrnpu_r2_rip_2_S2_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPK_r1_ --readFilesIn hnrnpk_hek293t_rip_1_S10_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPK_r2_ --readFilesIn hnrnpk_hek293t_rip_2_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_IGG_r1_ --readFilesIn igg_hek293t_rip_1_S13_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_IGG_r2_ --readFilesIn igg_hek293t_rip_2_S17_R1_001.fastq"

# merge replicates together for UCSC genome browser wiggle tracks

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r1_rip_5_S5_R1_001.fastq,G10_293T_flagtagged_GFP_RMCE_flag_r2_rip_6_S6_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r1_rip_5_S11_R1_001.fastq,G10_293T_flagtagged_GFP_RMCE_input_r2_rip_6_S12_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r1_rip_7_S7_R1_001.fastq,G10_293T_flagtagged_HNRNPK_RMCE_flag_r2_rip_8_S8_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r1_rip_7_S13_R1_001.fastq,G10_293T_flagtagged_HNRNPK_RMCE_input_r2_rip_8_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r1_rip_9_S9_R1_001.fastq,G10_293T_flagtagged_HNRNPU_RMCE_flag_r2_rip_10_S10_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r1_rip_9_S15_R1_001.fastq,G10_293T_flagtagged_HNRNPU_RMCE_input_r2_rip_10_S16_R1_001.fastq"



sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_flag_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_flag_r1_rip_5_S5_R1_001.fastq,G10_293T_flagtagged_GFP_RMCE_flag_r2_rip_6_S6_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_GFP_input_hg38_ --readFilesIn G10_293T_flagtagged_GFP_RMCE_input_r1_rip_5_S11_R1_001.fastq,G10_293T_flagtagged_GFP_RMCE_input_r2_rip_6_S12_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_flag_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_flag_r1_rip_7_S7_R1_001.fastq,G10_293T_flagtagged_HNRNPK_RMCE_flag_r2_rip_8_S8_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPK_input_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPK_RMCE_input_r1_rip_7_S13_R1_001.fastq,G10_293T_flagtagged_HNRNPK_RMCE_input_r2_rip_8_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_flag_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_flag_r1_rip_9_S9_R1_001.fastq,G10_293T_flagtagged_HNRNPU_RMCE_flag_r2_rip_10_S10_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix hm_HNRNPU_input_hg38_ --readFilesIn G10_293T_flagtagged_HNRNPU_RMCE_input_r1_rip_9_S15_R1_001.fastq,G10_293T_flagtagged_HNRNPU_RMCE_input_r2_rip_10_S16_R1_001.fastq"




sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPU_ --readFilesIn G10_293T_hnrnpu_r1_rip_1_S1_R1_001.fastq,G10_293T_hnrnpu_r2_rip_2_S2_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_HNRNPK_ --readFilesIn hnrnpk_hek293t_rip_1_S10_R1_001.fastq,hnrnpk_hek293t_rip_2_S14_R1_001.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix h_IGG_ --readFilesIn igg_hek293t_rip_1_S13_R1_001.fastq,igg_hek293t_rip_2_S17_R1_001.fastq"


### 1B. MACS peak-calling

mkdir bams
mv *_Aligned.out.sam ./bams/

cd ..
# split sam by strand
sbatch samtools_binarize_split_strands_1_30_24.sh unpaired reversed ./bams
# output saved under folder strands


# randomize strand
cd strands
### use merged files only, same as previous pipeline
sbatch run_rand.sh



# clean up names
for f in *_Aligned.out_neg_rand.bam; do
  [ -e "$f" ] || continue
  mv -i -- "$f" "${f/%_Aligned.out_neg_rand.bam/_neg_rand.bam}"
done

for f in *_Aligned.out_pos_rand.bam; do
  [ -e "$f" ] || continue
  mv -i -- "$f" "${f/%_Aligned.out_pos_rand.bam/_pos_rand.bam}"
done

# call peaks
sbatch run_macs.sh



#### Part 2: featureCounts count number of reads in each peak

## convert broadPeak to SAF file in preparation for featureCounts

for f in *_neg_peaks_ucsc.broadPeak; do
    echo "Processing $f ..."
    base="${f%_neg_peaks_ucsc.broadPeak}"
    cat "${base}_neg_peaks_ucsc.broadPeak" | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"-"}' > "${base}_saf.txt"
    cat "${base}_pos_peaks_ucsc.broadPeak" | sed 1d | awk -v OFS="\t" '{print $1,$2,$3,"+"}' >> "${base}_saf.txt"
    cat "${base}_saf.txt" | awk -v OFS="\t" '{print NR,$0}' | sed '1 i\GeneID\tChr\tStart\tEnd\tStrand' > "${base}_peaks.saf"
    
done



## featureCounts with SAF
# count reads in each peak for all samples with its proper control

cd ..
cd bams

sbatch run_fc.sh

