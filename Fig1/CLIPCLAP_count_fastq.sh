#!/bin/bash
#SBATCH -p general
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# Output file
out="CLIPCLAP_fastq_read_counts.txt"
echo -e "file\treads" > "$out"

# Make globs that don't match just disappear instead of staying literal
shopt -s nullglob

# Loop over fastq files
for f in *.fastq; do
    [ -e "$f" ] || continue  # skip if no match

    fname=$(basename "$f")

    lines=$(cat "$f" | wc -l)

    # 4 lines per read
    reads=$((lines / 4))

    echo -e "${fname}\t${reads}" >> "$out"
done
