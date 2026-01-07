# count hm samples with filtered human peaks

module load subread/2.1.1


featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPK_10xk.saf" -o "hm_HNRNPK_flag_merg_hg38_fc" hm_HNRNPK_flag_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPK_10xk.saf" -o "hm_HNRNPK_flag_r1_hg38_fc" hm_HNRNPK_flag_r1_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPK_10xk.saf" -o "hm_HNRNPK_flag_r2_hg38_fc" hm_HNRNPK_flag_r2_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPK_10xk.saf" -o "hm_HNRNPK_flag_ck_hg38_fc" hm_GFP_flag_hg38_Aligned.out.sam


featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPU_10xk.saf" -o "hm_HNRNPU_flag_merg_hg38_fc" hm_HNRNPU_flag_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPU_10xk.saf" -o "hm_HNRNPU_flag_r1_hg38_fc" hm_HNRNPU_flag_r1_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPU_10xk.saf" -o "hm_HNRNPU_flag_r2_hg38_fc" hm_HNRNPU_flag_r2_hg38_Aligned.out.sam

featureCounts -s 2 -F SAF -a "../saf/hg38_HNRNPU_10xk.saf" -o "hm_HNRNPU_flag_ck_hg38_fc" hm_GFP_flag_hg38_Aligned.out.sam



