import os, sys
import pandas as pd

contig_df = pd.read_csv("/Users/cmdb/Desktop/quant_assignment_1/asm/contigs.fasta.fai", sep = "\t", header  = None)
contig_length = contig_df.iloc[:, 1].to_list()

total_length = 0
for length in contig_length:
	total_length = total_length + int(length)

N50 = 0
sorted_contig_length = sorted(contig_length, reverse = True)

for lengths in sorted_contig_length:
	if N50 > total_length/2:
		break
	N50 = N50 + lengths

print(N50)