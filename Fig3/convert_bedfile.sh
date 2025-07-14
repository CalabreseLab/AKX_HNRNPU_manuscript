# write BED files based on $1_fc_rpm_2igg_ind.txt which identify the peaks with >2igg in >=2 reps
cat alyref_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > alyref_peaks_2igg_2reps.bed

cat ythdc1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > ythdc1_peaks_2igg_2reps.bed

cat hnrnpc_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > hnrnpc_peaks_2igg_2reps.bed

cat hnrnpk_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > hnrnpk_peaks_2igg_2reps.bed

cat hnrnpm_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > hnrnpm_peaks_2igg_2reps.bed

cat hnrnpu_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > hnrnpu_peaks_2igg_2reps.bed

cat hnrnpl_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > hnrnpl_peaks_2igg_2reps.bed

cat igg_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > igg_peaks_2igg_2reps.bed

cat input_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > input_peaks_2igg_2reps.bed

cat lbr_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > lbr_peaks_2igg_2reps.bed

cat matrin3_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > matrin3_peaks_2igg_2reps.bed

cat pabpn1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > pabpn1_peaks_2igg_2reps.bed

cat ptbp1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > ptbp1_peaks_2igg_2reps.bed

cat rbm15_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > rbm15_peaks_2igg_2reps.bed

cat rybp_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > rybp_peaks_2igg_2reps.bed

cat srsf1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > srsf1_peaks_2igg_2reps.bed

cat nudt21_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > nudt21_peaks_2igg_2reps.bed

cat supt16h_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > supt16h_peaks_2igg_2reps.bed

cat u2af35_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > u2af35_peaks_2igg_2reps.bed

cat u2af65_fc_rpm_2igg.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > u2af65_peaks_2igg_2reps.bed

cat xrn2_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > xrn2_peaks_2igg_2reps.bed

cat safb_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > safb_peaks_2igg_2reps.bed

cat spen_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > spen_peaks_2igg_2reps.bed

cat ring1b_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > ring1b_peaks_2igg_2reps.bed

cat tia1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > tia1_peaks_2igg_2reps.bed

cat ciz1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > ciz1_peaks_2igg_2reps.bed

cat khdrbs1_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > khdrbs1_peaks_2igg_2reps.bed

cat sfpq_fc_rpm_2igg_2reps.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > sfpq_peaks_2igg_2reps.bed

cat nxf1_fc_rpm_2igg.txt | awk -v OFS="\t" '{print $2,$3-1,$4,$1,0,$5}' > nxf1_peaks_2igg_2reps.bed

