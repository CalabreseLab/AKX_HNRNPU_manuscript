#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH -N 1
#SBATCH -n 8


module load ucsctools
module load deeptools

CHROMSIZES="mm10_chrNameLength.txt"


bigWigMerge diff_final_rbm15_peaks_pos.bw diff_final_safb_peaks_pos.bw diff_final_ythdc1_peaks_pos.bw diff_final_srsf1_peaks_pos.bw diff_final_nudt21_peaks_pos.bw diff_final_supt16h_peaks_pos.bw diff_final_u2af35_peaks_pos.bw diff_final_u2af65_peaks_pos.bw diff_final_khdrbs1_peaks_pos.bw diff_final_nxf1_peaks_pos.bw ../Airn_merged/Airn_27pro_c1_pos.bg

bigWigMerge diff_final_rbm15_peaks_neg.bw diff_final_safb_peaks_neg.bw diff_final_ythdc1_peaks_neg.bw diff_final_srsf1_peaks_neg.bw diff_final_nudt21_peaks_neg.bw diff_final_supt16h_peaks_neg.bw diff_final_u2af35_peaks_neg.bw diff_final_u2af65_peaks_neg.bw diff_final_khdrbs1_peaks_neg.bw diff_final_nxf1_peaks_neg.bw ../Airn_merged/Airn_27pro_c1_neg.bg

echo "merging completed c1"

bigWigMerge diff_final_hnrnpk_peaks_pos.bw diff_final_rybp_peaks_pos.bw diff_final_hnrnpu_peaks_pos.bw diff_final_pabpn1_peaks_pos.bw diff_final_ring1b_peaks_pos.bw diff_final_hnrnpl_peaks_pos.bw ../Airn_merged/Airn_27pro_c2_pos.bg

bigWigMerge diff_final_hnrnpk_peaks_neg.bw diff_final_rybp_peaks_neg.bw diff_final_hnrnpu_peaks_neg.bw diff_final_pabpn1_peaks_neg.bw diff_final_ring1b_peaks_neg.bw diff_final_hnrnpl_peaks_neg.bw ../Airn_merged/Airn_27pro_c2_neg.bg

echo "merging completed c2"

bigWigMerge diff_final_spen_peaks_pos.bw diff_final_hnrnpm_peaks_pos.bw diff_final_xrn2_peaks_pos.bw diff_final_tia1_peaks_pos.bw diff_final_sfpq_peaks_pos.bw ../Airn_merged/Airn_27pro_c3_pos.bg

bigWigMerge diff_final_spen_peaks_neg.bw diff_final_hnrnpm_peaks_neg.bw diff_final_xrn2_peaks_neg.bw diff_final_tia1_peaks_neg.bw diff_final_sfpq_peaks_neg.bw ../Airn_merged/Airn_27pro_c3_neg.bg

echo "merging completed c3"

bigWigMerge diff_final_ptbp1_peaks_pos.bw diff_final_lbr_peaks_pos.bw ../Airn_merged/Airn_27pro_c4_pos.bg

bigWigMerge diff_final_ptbp1_peaks_neg.bw diff_final_lbr_peaks_neg.bw ../Airn_merged/Airn_27pro_c4_neg.bg

echo "merging completed c4"


scp diff_final_alyref_peaks_*.bw ../Airn_merged/

mv ../Airn_merged/diff_final_alyref_peaks_pos.bw ../Airn_merged/Airn_27pro_c5_pos.bw

mv ../Airn_merged/diff_final_alyref_peaks_neg.bw ../Airn_merged/Airn_27pro_c5_neg.bw

scp diff_final_matrin3_peaks_*.bw ../Airn_merged/

mv ../Airn_merged/diff_final_matrin3_peaks_pos.bw ../Airn_merged/Airn_27pro_c6_pos.bw

mv ../Airn_merged/diff_final_matrin3_peaks_neg.bw ../Airn_merged/Airn_27pro_c6_neg.bw

scp diff_final_hnrnpc_peaks_*.bw ../Airn_merged/

mv ../Airn_merged/diff_final_hnrnpc_peaks_pos.bw ../Airn_merged/Airn_27pro_c7_pos.bw

mv ../Airn_merged/diff_final_hnrnpc_peaks_neg.bw ../Airn_merged/Airn_27pro_c7_neg.bw

scp diff_final_ciz1_peaks_*.bw ../Airn_merged/

mv ../Airn_merged/diff_final_ciz1_peaks_pos.bw ../Airn_merged/Airn_27pro_c8_pos.bw

mv ../Airn_merged/diff_final_ciz1_peaks_neg.bw ../Airn_merged/Airn_27pro_c8_neg.bw

cd ../Airn_merged

echo "start convert"

for f in Airn_27pro_*.bg; do
  sort -k1,1 -k2,2n "$f" > "${f%.bg}.sorted.bg" \
    && bedGraphToBigWig "${f%.bg}.sorted.bg" "$CHROMSIZES" "${f%.bg}.bw"
done

echo "finished convert"