#!/bin/bash
#SBATCH -p general
#SBATCH --mem=32G
#SBATCH --time=05:00:00


module load samtools/1.22


# Output file
out="CLIPCLAP_pair_counts.txt"
echo -e "sample\thuman_pairs\tmouse_pairs" > "$out"

for f in *_hg38_pairedonly.sam; do
    sample=$(basename "$f" .sam)

    # Properly paired (-f 2), exclude unmapped + secondary (-F 260)
    # Assumes SAM already has *only unique mappers* from STAR.

    mouse_pairs=$(samtools view -f 2 -F 260 "$f" \
        | awk '$3 ~ /_mm10$/ {c++} END{print int(c/2)}')

    human_pairs=$(samtools view -f 2 -F 260 "$f" \
        | awk '$3 !~ /_mm10$/ {c++} END{print int(c/2)}')

    echo -e "${sample}\t${human_pairs}\t${mouse_pairs}" >> "$out"
done


# the counts in human and mouse are the same if counting *_mm10_pairedonly.sam