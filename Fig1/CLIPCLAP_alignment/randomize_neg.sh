#!/bin/bash
#SBATCH -p general
#SBATCH --mem=64G
#SBATCH --time=05:00:00


module load samtools    

for f in *_neg.bam; do
    base="${f%_neg.bam}"
    echo "Randomizing $f -> ${base}_neg_rand.bam"
    samtools view -h "$f" \
      | perl randomize_pairs.pl \
      | samtools view -b -o "${base}_neg_rand.bam" \
      2> "${base}_neg_rand.log"
done