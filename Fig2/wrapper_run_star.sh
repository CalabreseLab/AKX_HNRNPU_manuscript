#!/bin/bash
#SBATCH --output=star.out
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3:00:00
#SBATCH --mem=32g

module load star

for faq in *.fastq
do
  base=${faq%.fastq}
  sbatch --export=ALL,base=${base} run_star.sh ${base}
done
