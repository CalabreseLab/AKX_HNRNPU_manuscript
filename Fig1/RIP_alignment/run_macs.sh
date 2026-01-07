#!/bin/bash
#SBATCH -p general
#SBATCH --mem=96G
#SBATCH --time=2-0


module load macs/2.2.9.1

# -n: The prefix string for output files

for f in *_rand.bam; do
    base="${f%_rand.bam}"
    echo "calling peaks $f"
    macs2 callpeak -t "${base}_rand.bam" --broad --broad-cutoff 0.3 --max-gap 76 --outdir "${base}_peaks" -n "${base}"
    # add trackline for UCSC
    cat "${base}_peaks/${base}_peaks.broadPeak" | sed "1 i\track type=broadPeak name="${base}_peaks" description="${base}_peaks" nextItemButton=on" > "${base}_peaks_ucsc.broadPeak"


done


