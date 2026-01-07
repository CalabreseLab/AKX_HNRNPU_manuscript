#!/bin/bash
#SBATCH -p general
#SBATCH --mem=64G
#SBATCH --time=05:00:00


module load samtools    

for f in *_pos.bam; do
    base="${f%_pos.bam}"
    echo "Randomizing $f -> ${base}_pos_rand.bam"
    samtools view -h "$f" \
      | perl randomize_pairs.pl \
      | samtools view -b -o "${base}_pos_rand.bam" \
      2> "${base}_pos_rand.log"
done
