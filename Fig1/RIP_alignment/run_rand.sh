#!/bin/bash
#SBATCH -p general
#SBATCH --mem=96G
#SBATCH --time=1-0



for f in !(*_r1*|*_r2*|*_input*).bam; do
    base="${f%.bam}"
    echo "Randomizing $f -> ${base}_rand.bam"

    samtools view -h "$f" \
      | perl macs_strand_rand_sam_12082025.pl \
      | samtools view -b -o "${base}_rand.bam" \
      2> "${base}_rand.log"

done

