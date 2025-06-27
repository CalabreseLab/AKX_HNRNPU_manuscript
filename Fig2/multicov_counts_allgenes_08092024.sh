#!/bin/bash

module load bedtools

# exclude input and other 6 rbps from the analysis
ls *_sorted.bam | grep -Ev 'input|bmi1|epop|ezh2|jarid2|suz12|mtf2' > multicov_files_list_allgenes_08092024.txt

#ls *_sorted.bam | grep -v 'input' > multicov_files_list_allgenes.txt
bam_files=$(cat multicov_files_list_allgenes_08092024.txt | tr '\n' ' ')

# ls *_sorted.bam > multicov_files_list_allgenes.txt
# bam_files=$(cat multicov_files_list_allgenes.txt | tr '\n' ' ')

for file in ../bedfiles/*
do
    base=$(basename "$file" .bed)
    sbatch -p general --time=1-0 -n 8 -N 1 --mem=64g -o "/work/users/s/h/shuang9/rip/multicov_out/${base}_multicov.out" --wrap="bedtools multicov -S -D -q 30 -bams $bam_files -bed $file"
done
