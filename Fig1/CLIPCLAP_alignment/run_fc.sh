#!/bin/bash
#SBATCH -p general
#SBATCH --mem=96G
#SBATCH --time=1-0


module load subread/2.1.1

# count human
# when we call peaks in human, the - tag will be our control for +tag.

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLAP_hg38_fc" SAFA_PlusTag_CLAP_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLAP_hg38_ck_fc" SAFA_MinusTag_CLAP_hg38_pairedonly.sam


featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLAP_Rep1_hg38_fc" SAFA_PlusTag_CLAP_Rep1_hg38_pairedonly.sam


featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLAP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLAP_Rep2_hg38_fc" SAFA_PlusTag_CLAP_Rep2_hg38_pairedonly.sam




featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLIP_hg38_fc" SAFA_PlusTag_CLIP_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLIP_hg38_ck_fc" SAFA_MinusTag_CLIP_hg38_pairedonly.sam


featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLIP_Rep1_hg38_fc" SAFA_PlusTag_CLIP_Rep1_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_PlusTag_CLIP_hg38_pairedonly_peaks.saf" -o "SAFA_PlusTag_CLIP_Rep2_hg38_fc" SAFA_PlusTag_CLIP_Rep2_hg38_pairedonly.sam





# count mouse
# when we call peaks in mouse, the +tag will be the control for -tag

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_MinusTag_CLAP_mm10_pairedonly_peaks.saf" -o "SAFA_MinusTag_CLAP_mm10_fc" SAFA_MinusTag_CLAP_mm10_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_MinusTag_CLAP_mm10_pairedonly_peaks.saf" -o "SAFA_MinusTag_CLAP_mm10_ck_fc" SAFA_PlusTag_CLAP_mm10_pairedonly.sam


featureCounts -s 2 -p -F SAF -a "../strands/SAFA_MinusTag_CLIP_mm10_pairedonly_peaks.saf" -o "SAFA_MinusTag_CLIP_mm10_fc" SAFA_MinusTag_CLIP_mm10_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../strands/SAFA_MinusTag_CLIP_mm10_pairedonly_peaks.saf" -o "SAFA_MinusTag_CLIP_mm10_ck_fc" SAFA_PlusTag_CLIP_mm10_pairedonly.sam
