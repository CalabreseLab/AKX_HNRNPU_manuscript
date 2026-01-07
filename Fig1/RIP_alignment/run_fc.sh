#!/bin/bash
#SBATCH -p general
#SBATCH --mem=96G
#SBATCH --time=1-0


module load subread/2.1.1

# count human only cells
# IgG is the control

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPK_peaks.saf" -o "h_HNRNPK_r1_fc" h_HNRNPK_r1_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPK_peaks.saf" -o "h_HNRNPK_r2_fc" h_HNRNPK_r2_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPK_peaks.saf" -o "h_HNRNPK_merg_fc" h_HNRNPK_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPK_peaks.saf" -o "h_HNRNPK_ck_fc" h_IGG_Aligned.out.sam


featureCounts -s 2 -F SAF -a "../strands/h_HNRNPU_peaks.saf" -o "h_HNRNPU_r1_fc" h_HNRNPU_r1_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPU_peaks.saf" -o "h_HNRNPU_r2_fc" h_HNRNPU_r2_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPU_peaks.saf" -o "h_HNRNPU_merg_fc" h_HNRNPU_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/h_HNRNPU_peaks.saf" -o "h_HNRNPU_ck_fc" h_IGG_Aligned.out.sam



# count hm samples in mouse use peaks defined by hm samples
# the FLAG-GFP will be the control
featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPK_flag_peaks.saf" -o "hm_HNRNPK_flag_merg_fc" hm_HNRNPK_flag_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPK_flag_peaks.saf" -o "hm_HNRNPK_flag_r1_fc" hm_HNRNPK_flag_r1_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPK_flag_peaks.saf" -o "hm_HNRNPK_flag_r2_fc" hm_HNRNPK_flag_r2_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPK_flag_peaks.saf" -o "hm_HNRNPK_flag_ck_fc" hm_GFP_flag_Aligned.out.sam


featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPU_flag_peaks.saf" -o "hm_HNRNPU_flag_merg_fc" hm_HNRNPU_flag_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPU_flag_peaks.saf" -o "hm_HNRNPU_flag_r1_fc" hm_HNRNPU_flag_r1_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPU_flag_peaks.saf" -o "hm_HNRNPU_flag_r2_fc" hm_HNRNPU_flag_r2_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../strands/hm_HNRNPU_flag_peaks.saf" -o "hm_HNRNPU_flag_ck_fc" hm_GFP_flag_Aligned.out.sam


# hm samples in human use peaks defined by h_ samples (human only) will be counted in a later step
# after filtering the human only peaks




