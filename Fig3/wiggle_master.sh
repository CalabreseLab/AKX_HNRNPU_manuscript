# master pipeline for generating wiggle tracks

module load samtools
module load bedtools


bash bam_to_bed12_1_30_24.sh ./bams


bash run_count_reads_1_30_24.sh ./bams

# copy over mm10 chrName and chrLength file 
paste -d $'\t' chrName.txt chrLength.txt > mm10_chrNameLength.txt

# add up the neg and pos counts for each RIP and IgG-RIP and save as *_peaks_readcounts.txt
cd readcounts
bash sum_readcounts.sh
rm *_neg_readcounts.txt
rm *_pos_readcounts.txt

cd ..
sbatch make_wiggle_script_input_10_24_24.sh bedfiles mm10_chrNameLength.txt stranded n 50 readcounts/ n

bash run_wiggle_script.sh

mv *.wig ./wigs

module load ucsctools

bash make_bigwigs_1_30_24.sh mm10_chrNameLength.txt ./wigs

# Subtract normalized IgG signal from RIP signal in wiggle, to generate bkg corrected wiggle
module load ucsctools
module load deeptools/3.5.1

cd bigWigoutputs
sbatch subtract_bigwig.sh


