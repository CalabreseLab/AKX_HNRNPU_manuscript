

module load bedtools

for f in *.bed; do bedtools merge -s -d 100 -c 4,5,6 -o distinct,sum,distinct -i "$f" > "${f%.bed}_merged100.bed"; done

for f in *_merged100.bed; do bedtools getfasta -s -name -fi /mouse/mm10_genome_filtered.fa -bed "$f" -fo "${f%.bed}.fa"; done


# concatenate all merged peaks into one bedfile and sort for XAK

cat Xist_*_merged100.bed | sort -k1,1 -k2,2n > Xist_all_merged100.sorted.bed

cat Airn_*_merged100.bed | sort -k1,1 -k2,2n > Airn_all_merged100.sorted.bed

cat Kot1_*_merged100.bed | sort -k1,1 -k2,2n > Kot1_all_merged100.sorted.bed


bedtools subtract -s -a Xist.bed -b Xist_all_merged100.sorted.bed > "Xist_null.bed"

bedtools subtract -s -a Airn.bed -b Airn_all_merged100.sorted.bed > "Airn_null.bed"

bedtools subtract -s -a Kot1.bed -b Kot1_all_merged100.sorted.bed > "Kot1_null.bed"



bedtools getfasta -s -name -fi /mouse/mm10_genome_filtered.fa -bed Xist_null.bed -fo "Xist_null.fa"

bedtools getfasta -s -name -fi /mouse/mm10_genome_filtered.fa -bed Airn_null.bed -fo "Airn_null.fa"

bedtools getfasta -s -name -fi /mouse/mm10_genome_filtered.fa -bed Kot1_null.bed -fo "Kot1_null.fa"





