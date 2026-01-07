#!/bin/bash
#SBATCH -p general
#SBATCH --mem=96G
#SBATCH --time=2-0


module load macs/2.2.9.1

# -n: The prefix string for output files
# although it is recommended: macs2 peak calling with paired end reads -f BAMPE
# For CLAP-style data, it’s usually better to stay in the “BAM” (single-end) mode (default or -f BAM) using the 5' ends rather than BAMPE.
# with no -f, MACS2 treated your BAM as single-end “BAM”, even though the file contains paired-end reads
# can see this in the log output
# In this mode MACS2:
# Uses each read’s 5' end as the tag position.
# Builds a shift/extension model from the plus/minus strand cross-correlation:
# It looked for “paired plus/minus strand peaks…”
# Estimated a fragment length d = 84 bp.
# It then extends each 5' tag to length d to approximate the real fragment and calls peaks from that coverage.
# your wiggle track is almost certainly based on raw read coverage or 5' positions, not on reconstructed paired-end fragments, so the MACS peaks from this mode line up nicely with what you see.

# In BAMPE mode, MACS2 no longer:infers d from cross-correlation in the same way, or uses only the 5' end of each read.
# Instead, for each properly paired fragment: It finds the template span between mate1 and mate2 (the insert).
# Treats that entire fragment as a single tag (one fragment = one observation).
# Coverage is constructed from fragment spans, not individual read ends.

# So now: Your peak calls reflect fragment spans and/or fragment centers, not 5′ read piles.
# If your wiggle track is built from: just one mate, or all reads treated independently (ignoring pairing),
# it will not overlap well with fragment-based peaks. Often the entire signal looks shifted, broadened, or “moved” relative to the raw read track.
# For CLAP / CLIP-like assays, this mismatch can be dramatic, because biologically meaningful positions are usually at the crosslink 5′ ends, not over the full cDNA fragment.

# For CLIP/CLAP/RBP-binding experiments:The 5′ end marks the crosslink site.
# The fragment insert can be quite variable and may extend far from the crosslink.
# Using BAMPE causes MACS2 to use the full fragment spans, smearing signal across broader regions and shifting the apparent “peak” away from the sharp pileup of 5′ ends.
# So, BAM mode aligns with biology here; BAMPE does not.



for f in *_rand.bam; do
    base="${f%_rand.bam}"
    echo "calling peaks $f"
    macs2 callpeak -t "${base}_rand.bam" --broad --broad-cutoff 0.3 --max-gap 76 --outdir "${base}_peaks" -n "${base}"
    # add trackline for UCSC
    cat "${base}_peaks/${base}_peaks.broadPeak" | sed "1 i\track type=broadPeak name="${base}_peaks" description="${base}_peaks" nextItemButton=on" > "${base}_peaks_ucsc.broadPeak"


done


