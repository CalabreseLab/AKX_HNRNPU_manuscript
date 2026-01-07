# count CLIP with CLAP bedfile

module load subread/2.1.1

featureCounts -s 2 -p -F SAF -a "../saf/hg38_CLAP_10xk.saf" -o "hg38CLAP10k_CLIP_fc" SAFA_PlusTag_CLIP_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/hg38_CLAP_10xk.saf" -o "hg38CLAP10k_CLIP_ck_fc" SAFA_MinusTag_CLIP_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/hg38_CLAP_10xk.saf" -o "hg38CLAP10k_CLIP_Rep1_fc" SAFA_PlusTag_CLIP_Rep1_hg38_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/hg38_CLAP_10xk.saf" -o "hg38CLAP10k_CLIP_Rep2_fc" SAFA_PlusTag_CLIP_Rep2_hg38_pairedonly.sam



featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLIP_fc" SAFA_MinusTag_CLIP_mm10_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLIP_ck_fc" SAFA_PlusTag_CLIP_mm10_pairedonly.sam



##### check the replciates counts
featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLIP_ck1_fc" SAFA_PlusTag_CLIP_Rep1_mm10_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLIP_ck2_fc" SAFA_PlusTag_CLIP_Rep2_mm10_pairedonly.sam



featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLAP_ck1_fc" SAFA_PlusTag_CLAP_Rep1_mm10_pairedonly.sam

featureCounts -s 2 -p -F SAF -a "../saf/mm10_CLAP_10xk.saf" -o "mm10CLAP10k_CLAP_ck2_fc" SAFA_PlusTag_CLAP_Rep2_mm10_pairedonly.sam




