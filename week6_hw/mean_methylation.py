import sys, os
import pandas as pd

def get_dataframe(filename):
	df = pd.read_csv(filename, sep='\t', header = None)
	return(df)

#get a dictionary of gene names (key) and gene length (values are [start, end])
gene_name_and_location_dic = {}
bed_ref = get_dataframe("mm10_refseq_genes_chr6_50M_60M.bed")

col_size = bed_ref[12].size

for i in range(col_size):
	if bed_ref[12][i] in gene_name_and_location_dic:
		continue
	else:
		gene_name_and_location_dic[bed_ref[12][i]] = [bed_ref[4][i], bed_ref[5][i]]

#calculate mean methylation for each gene for each E4/E5.5 bedGraph result
e4 = pd.read_csv("SRR3083926_1.chr6_bismark_bt2_pe.sorted.bedGraph", sep='\t', header = None, skiprows = [0])
e5_5 = pd.read_csv("SRR3083929_1.chr6_bismark_bt2_pe.sorted.bedGraph", sep='\t', header = None, skiprows = [0])

def calc_mean_methylation(df, gene_dic):
	mean_methylation_dic = {}
	for gene in gene_name_and_location_dic:
		start = gene_name_and_location_dic[gene][0]
		end = gene_name_and_location_dic[gene][1]
		gene_segment = df[((df[1] >= start) & (df[2] <= end))]
		mean_methylation = gene_segment[3].divide(100).sum()
		mean_methylation_dic[gene] = mean_methylation
	return(mean_methylation_dic)

e4_mean_dic = calc_mean_methylation(e4, gene_name_and_location_dic)
e5_mean_dic = calc_mean_methylation(e5_5, gene_name_and_location_dic)

fold_change_dic = {}
for gene in e4_mean_dic:
	if e4_mean_dic[gene] == 0:
		continue
	if e4_mean_dic[gene] < e5_mean_dic[gene]:
		fold = e5_mean_dic[gene]/e4_mean_dic[gene]
		fold_change_dic[gene] = fold
	if e5_mean_dic[gene] < e4_mean_dic[gene]:
		fold = (e5_mean_dic[gene]/e4_mean_dic[gene]) * -1
		fold_change_dic[gene] = fold

with open("fold_change_e4_to_e5_5.txt", "a+") as newFile:
	for gene in fold_change_dic:
		newFile.write(gene + "\t" + str(fold_change_dic[gene]) + "\n")
