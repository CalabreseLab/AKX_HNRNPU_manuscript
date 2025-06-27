#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBACTH --mem=32g

module load samtools

for sam_file in *.sam
do
  base=${sam_file%_Aligned.out.sam}
  sbatch --export=ALL,base=${base} convert_sam_to_bam.sh ${base}
done
