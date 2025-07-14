#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8


module load ucsctools
module load deeptools

CHROMSIZES="mm10_chrNameLength.txt"
sample="$1"
suffix="$2"
base=${sample%_${suffix}.bw}
ctrl="igg_${base}_${suffix}.bw"

if [[ ! -e "$ctrl" ]]; then
  echo "Missing control for $sample: $ctrl not found" >&2
  continue
fi

echo "Processing $sample  (ctrl: $ctrl)"
bigwigCompare \
  --operation subtract \
  -b1 "$sample" \
  -b2 "$ctrl" \
  -o raw_diff_${base}_${suffix}.bw


echo "subtract completed"

bigWigToBedGraph raw_diff_${base}_${suffix}.bw raw_diff_${base}_${suffix}.bg
awk '{ if($4<0) $4=0; print }' raw_diff_${base}_${suffix}.bg > diff_pos_${base}_${suffix}.bg

echo "convert negtives completed"

bedGraphToBigWig diff_pos_${base}_${suffix}.bg "$CHROMSIZES" diff_final_${base}_${suffix}.bw

echo "Wrote diff_final_${base}_${suffix}.bw"
