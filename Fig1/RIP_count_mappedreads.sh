#!/bin/bash
#SBATCH -p general
#SBATCH --mem=32G
#SBATCH --time=05:00:00


module load samtools/1.22


# Output file
out="RIP_unimap_counts.txt"
echo -e "sample\thuman_pairs\tmouse_pairs" > "$out"

for f in bams/hm_*.sam; do
    sample=$(basename "$f" .sam)

    # Properly paired (-f 2), exclude unmapped + secondary (-F 260)
    # Assumes SAM already has *only unique mappers* from STAR.

    mouse_pairs=$(samtools view -F 4 "$f" \
        | awk '$3 !~ /_hg38$/ {c++} END{print c}')

    human_pairs=$(samtools view -F 4 "$f" \
        | awk '$3 ~ /_hg38$/ {c++} END{print c}')

    echo -e "${sample}\t${human_pairs}\t${mouse_pairs}" >> "$out"
done

