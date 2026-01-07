#!/bin/bash

#SBATCH -p general
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-0


### for all RIP seq samples, use STAR featureCounts to align and count to ERCC sequences 
### calculate RPM using total sequence count from fastq

## load modules
module load star/2.7.7a
module load samtools/1.22
module load subread/2.1.1
module load bedtools/2.30.0


#### Part 1: from rip, STAR genome alignment 
#### Part 2: featureCounts count number of reads 
#### Part 3: Convert to rpm = (count * 1000,000 / reads)


mkdir $1_ercc
cd $1_ercc

### 1. STAR alignment

## STAR alignment of all rips
# move rip files here
mkdir $1
cp ${2} ${1}
gunzip $1/*.fastq.gz

# use STAR aligner to align reads to ERCC.
rip_filenames=$(ls -m $1/*.fastq | sed -r 's/\s+//g' | tr -d '\n')

## STAR alignment of individual rips
# decide if there are more than 1 rip
rip_file_num=$(ls $1 | wc -l)
if [ $rip_file_num -gt 1 ]
then

    IFS=',' read -r -a rip_arr <<< "$rip_filenames" # make array with rip file names

    # align each rip file
    for i in ${!rip_arr[@]}; do # loop on indices of array
        rip_name="${rip_arr[$i]}" # for each file name
        echo The rip${i} file is $rip_name. >> $1_rip_info.txt # NOTE: check in $1_rip_info.txt to see which rip gets assigned to which number!
        j=$(($i + 1)) # index+1 so that it's 1 based for the file name

        # align to mm10 genome
        STAR --genomeDir /rip/rip_ercc/ERCC_STAR_index --runThreadN 8 --outFilterMultimapNmax 1 --outFileNamePrefix $1_rip${j}_ --readFilesIn $rip_name
        

        # count read for each individual rip
        featureCounts -s 2 -g gene_id -a /rip/rip_ercc/ERCC92.gtf -o $1_rip${j}_fc $1_rip${j}_Aligned.out.sam

        # extract raw reads

        sed 1,2d $1_rip${j}_fc | cut -f7 > $1_rip${j}_counts.txt

        # rpm conversion for individual rips:
        # count reads
        fastq_lin_num=$(wc -l ${rip_arr[$i]} | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
        (( reads = fastq_lin_num / 4 ))
        # rpm conversion
        cat $1_rip${j}_counts.txt | awk -v reads=$reads '{print (($1*1000000)/reads)}' > $1_rip${j}_rpm.txt

        # prep for concatenate rpm files and file header
        if [ $j == 1 ]
        then
            rip_count_header="$1_rip1_counts"
            rip_rpm_header="$1_rip1_rpm"
            rip_count_filename="$1_rip1_counts.txt"
            rip_rpm_filename="$1_rip1_rpm.txt"
        else
            rip_count_header="${rip_count_header}\t$1_rip${j}_counts"
            rip_rpm_header="${rip_rpm_header}\t$1_rip${j}_rpm"
            rip_count_filename="${rip_count_filename} $1_rip${j}_counts.txt"
            rip_rpm_filename="${rip_rpm_filename} $1_rip${j}_rpm.txt"
        fi

    done
    # extract the first 6 col of fc to use as the features
    sed 1,2d $1_rip1_fc | cut -f1-6 > $1_features.txt
    # make summary of rpm
    paste $1_features.txt $rip_count_filename $rip_rpm_filename | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t${rip_count_header}\t${rip_rpm_header}" > $1_fc_rpm.txt


else
    STAR --genomeDir /rip/rip_ercc/ERCC_STAR_index --runThreadN 8 --outFilterMultimapNmax 1 --outFileNamePrefix $1_ --readFilesIn $rip_filenames
    #samtools view -h -Sq 30 $1_Aligned.out.sam > $1_Aligned_filteredsq30.out.sam
    featureCounts -s 2 -g gene_id -a /rip/rip_ercc/ERCC92.gtf -o $1_fc $1_Aligned.out.sam
    # extract raw reads
    sed 1,2d $1_fc | cut -f7 > $1_counts.txt
    fastq_lin_num=$(wc -l $1/*$1*.fastq | tail -1 | sed 's/^ *//g' | cut -d ' ' -f1)
    (( reads = fastq_lin_num / 4 ))
    # rpm conversion
    cat $1_counts.txt | awk -v reads=$reads '{print (($1*1000000)/reads)}' > $1_rpm.txt
    sed 1,2d $1_fc | cut -f1-6 > $1_features.txt
    paste $1_features.txt $1_counts.txt $1_rpm.txt | sed "1 i\GeneID\tChr\tStart\tEnd\tStrand\tLength\t$1_counts\t$1_rpm" > $1_fc_rpm.txt

fi


