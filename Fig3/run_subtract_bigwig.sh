#!/bin/bash
shopt -s nullglob

for suffix in neg pos; do
  for sample in *_peaks_${suffix}.bw; do
    # skip pure‐control files
    [[ $(basename "$sample") == igg_* ]] && continue

    # submit one job per sample+suffix
    sbatch subtract_bigwig.sh "$sample" "$suffix"
  done
done
