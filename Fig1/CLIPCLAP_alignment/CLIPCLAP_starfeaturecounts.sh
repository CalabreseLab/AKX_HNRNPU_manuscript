#!/bin/bash

# when we call peaks in human, the (-) will be our control for (+). 
# But when we call peaks in mouse, the (+) will be the control for (-)

## load modules
module load star//2.7.11b
module load samtools/1.22
module load macs/2.2.9.1
module load subread/2.1.1
module load bedtools/2.31.1


#### Part 1: from rip, STAR genome alignment and MACS peak-calling 

### 1A. STAR alignment

## STAR alignment of all samples

# use STAR aligner to align reads to human genome
sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_CLAP_hg38_ --readFilesIn SAFA_MinusTag_CLAP_r1.fastq SAFA_MinusTag_CLAP_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_CLIP_hg38_ --readFilesIn SAFA_MinusTag_CLIP_r1.fastq SAFA_MinusTag_CLIP_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_Input_hg38_ --readFilesIn SAFA_MinusTag_Input_r1.fastq SAFA_MinusTag_Input_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_hg38_ --readFilesIn SAFA_PlusTag_CLAP_Rep1_r1.fastq,SAFA_PlusTag_CLAP_Rep2_r1.fastq SAFA_PlusTag_CLAP_Rep1_r2.fastq,SAFA_PlusTag_CLAP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_hg38_ --readFilesIn SAFA_PlusTag_CLIP_Rep1_r1.fastq,SAFA_PlusTag_CLIP_Rep2_r1.fastq SAFA_PlusTag_CLIP_Rep1_r2.fastq,SAFA_PlusTag_CLIP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_Input_hg38_ --readFilesIn SAFA_PlusTag_Input_r1.fastq SAFA_PlusTag_Input_r2.fastq"


# use STAR aligner to align reads to mouse genome

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_CLAP_mm10_ --readFilesIn SAFA_MinusTag_CLAP_r1.fastq SAFA_MinusTag_CLAP_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_CLIP_mm10_ --readFilesIn SAFA_MinusTag_CLIP_r1.fastq SAFA_MinusTag_CLIP_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_MinusTag_Input_mm10_ --readFilesIn SAFA_MinusTag_Input_r1.fastq SAFA_MinusTag_Input_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_mm10_ --readFilesIn SAFA_PlusTag_CLAP_Rep1_r1.fastq,SAFA_PlusTag_CLAP_Rep2_r1.fastq SAFA_PlusTag_CLAP_Rep1_r2.fastq,SAFA_PlusTag_CLAP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_mm10_ --readFilesIn SAFA_PlusTag_CLIP_Rep1_r1.fastq,SAFA_PlusTag_CLIP_Rep2_r1.fastq SAFA_PlusTag_CLIP_Rep1_r2.fastq,SAFA_PlusTag_CLIP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_Input_mm10_ --readFilesIn SAFA_PlusTag_Input_r1.fastq SAFA_PlusTag_Input_r2.fastq"


## STAR alignment of individual rips

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_Rep1_hg38_ --readFilesIn SAFA_PlusTag_CLAP_Rep1_r1.fastq SAFA_PlusTag_CLAP_Rep1_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_Rep2_hg38_ --readFilesIn SAFA_PlusTag_CLAP_Rep2_r1.fastq SAFA_PlusTag_CLAP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_Rep1_hg38_ --readFilesIn SAFA_PlusTag_CLIP_Rep1_r1.fastq SAFA_PlusTag_CLIP_Rep1_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_hg38_mm10_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_Rep2_hg38_ --readFilesIn SAFA_PlusTag_CLIP_Rep2_r1.fastq SAFA_PlusTag_CLIP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_Rep1_mm10_ --readFilesIn SAFA_PlusTag_CLAP_Rep1_r1.fastq SAFA_PlusTag_CLAP_Rep1_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLAP_Rep2_mm10_ --readFilesIn SAFA_PlusTag_CLAP_Rep2_r1.fastq SAFA_PlusTag_CLAP_Rep2_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_Rep1_mm10_ --readFilesIn SAFA_PlusTag_CLIP_Rep1_r1.fastq SAFA_PlusTag_CLIP_Rep1_r2.fastq"

sbatch -p general --time=10:00:00 -n 8 -N 1 --mem=96g --wrap="STAR --genomeDir STAR_mm10_hg38_index --runThreadN 8 --outSAMtype SAM --outFilterMultimapNmax 1 --outFileNamePrefix SAFA_PlusTag_CLIP_Rep2_mm10_ --readFilesIn SAFA_PlusTag_CLIP_Rep2_r1.fastq SAFA_PlusTag_CLIP_Rep2_r2.fastq"


######
# keep only properly paired reads

for f in *_Aligned.out.sam; do
    echo "Processing $f ..."
    samtools view -h -f 2 "$f" > "${f%_Aligned.out.sam}_pairedonly.sam"
done

### 1B. MACS peak-calling

mkdir bams
mv *_pairedonly.sam ./bams/

cd ..
# split sam by strand
sbatch samtools_binarize_split_strands_1_30_24.sh paired reversed ./bams
# output saved under folder strands


# randomize strand
cd strands

# submit job using:
sbatch randomize_pos.sh
sbatch randomize_neg.sh
# output is randomized bam file _pos/neg_rand.bam


# although for paired end seq it is recommended to use -f BAMPE
# For CLAP-style data, it’s usually better to stay in the “BAM” (single-end) mode (default or -f BAM) using the 5' ends rather than BAMPE.
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
