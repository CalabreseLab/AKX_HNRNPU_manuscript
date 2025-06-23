# Written on 2/27/23 to calculate the nucleotide content of a given genome
# and to develop weighted randomized control sequences the same length as
# given query sequences (originally for Rachel SAFB motifs)

import numpy as np
import os

A_count = 0
T_count = 0
G_count = 0
C_count = 0
other_count = 0

with open("/work/users/s/h/shuang9/rip/GRCm38.primary_assembly.genome.fa", "r") as genomefile:
    for line in genomefile:
        if line[0] != ">":
            for char in line:
                if char == "A":
                    A_count += 1
                elif char == "T":
                    T_count += 1
                elif char == "G":
                    G_count += 1
                elif char == "C":
                    C_count += 1
                    
print(A_count)
print(T_count)
print(G_count)
print(C_count)

total_nts = A_count + T_count + G_count + C_count
A_frac = A_count/total_nts
T_frac = T_count/total_nts
G_frac = G_count/total_nts
C_frac = C_count/total_nts

nts = ["A","T","G","C"]

def create_background_seq_file(controlfile, peakseqfile):
    counter = 0
    for line in peakseqfile:
        if line[0] != ">":
            counter += 1
            outcontrolfile.write(">%s\n" % counter)
            controlfile.write(str("%s\n" % "".join(np.random.choice(nts,len(line.strip()),p=[A_frac,T_frac,G_frac,C_frac]))))


# get all the fasta files under top1kfasta folder
fas = os.listdir("/work/users/s/h/shuang9/rip/top1kfasta/")
fas = [entry for entry in fas if entry.endswith('.fa')]

for fa in fas:
    with open("/work/users/s/h/shuang9/rip/top1kfasta/%s" % fa, "r") as inpeakfile:
        with open("/work/users/s/h/shuang9/rip/top1kcontrol/control_%s" % fa, "w") as outcontrolfile:
            create_background_seq_file(outcontrolfile, inpeakfile)



