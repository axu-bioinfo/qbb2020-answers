import os, sys
import pandas as pd

#Run this python script by changing the file path for file BDGP6.Ensembl.81.gtf (not the file name, the actual path)
#in the last two lines
#This code is ugly so I'll add the output to the end of the file as a comment so you know it actually
#works

def read_gtf_df(file_path, file_name):
	just_gene_data = os.system("grep -v \'#\' " + file_path + file_name + " | grep \'3R\t\' > " + file_path + "temp.txt")
	final_3R_df = pd.read_csv(file_path + "temp.txt", sep = "\t", header = None)
	os.system("rm " + file_path + "temp.txt")
	return(final_3R_df)

def binary_search(low, mid, high, df, pos_num):
	start = df.iloc[mid][3]
	end = df.iloc[mid][4]
	
	if pos_num > start and pos_num < end:
		return("found it")

	if pos_num < start:
		return("lower")

	if pos_num > end:
		return("higher")


def run_binary_search(position_num, gtf_df):
	all_3R_prot = len(gtf_df.index)
	all_3R_prot_index = all_3R_prot - 1

	next_3R_prot_index = all_3R_prot_index//2
	low = 0
	high = all_3R_prot_index
	final_index = 0
	count_iter = 0
	while True:
		result = binary_search(low, next_3R_prot_index, high, gtf_df, position_num)
		if result == "found it":
			final_index = next_3R_prot_index
			break
		if result == "lower":
			if next_3R_prot_index == high:
				last_end = gtf_df.iloc[next_3R_prot_index][4]
				next_begin = gtf_df.iloc[next_3R_prot_index + 1][3]
				if abs(position_num - last_end) < abs(position_num - next_begin):
					final_index = next_3R_prot_index
				if abs(position_num - last_end) > abs(position_num - next_begin):
					final_index = next_3R_prot_index + 1
				break
			high = next_3R_prot_index
			next_3R_prot_index = (low + next_3R_prot_index)//2
		if result == "higher":
			if next_3R_prot_index == low:
				last_end = gtf_df.iloc[next_3R_prot_index][4]
				next_begin = gtf_df.iloc[next_3R_prot_index + 1][3]
				if abs(position_num - last_end) < abs(position_num - next_begin):
					final_index = next_3R_prot_index
				if abs(position_num - last_end) > abs(position_num - next_begin):
					final_index = next_3R_prot_index + 1
				break
			low = next_3R_prot_index
			next_3R_prot_index = (high + next_3R_prot_index)//2
		count_iter = count_iter + 1
	print("This thing iterated " + str(count_iter) + " times before finding the answer")
	return(final_index)

def get_the_output(num, df, close_index):
	gene_name_start = df.iloc[close_index][8].find("gene_name \"") + 11
	gene_name_end = df.iloc[close_index][8].find("\"", gene_name_start)
	gene_name = df.iloc[close_index][8][gene_name_start:gene_name_end]
	print("The gene name is " + gene_name)

	gene_id_start = df.iloc[close_index][8].find("gene_id \"") + 9
	gene_id_end = df.iloc[close_index][8].find("\"", gene_id_start)
	gene_id = df.iloc[close_index][8][gene_id_start:gene_id_end]
	print("The gene id is " + gene_id)

	biotype_start = df.iloc[close_index][8].find("gene_biotype \"") + 14
	biotype_end = df.iloc[close_index][8].find("\"", biotype_start)
	biotype = df.iloc[close_index][8][biotype_start:biotype_end]
	print("The biotype is " + biotype)

	begin = df.iloc[close_index][3]
	end = df.iloc[close_index][4]

	if abs(begin - num) < abs(end - num):
		print("Shortest distance from the beginning of gene is " + str(abs(begin - num)) + " bp")
	else:
		print("Shortest distance from the end of gene is " + str(abs(end - num)) + " bp")




closest_index = run_binary_search(21378950, read_gtf_df("/Users/cmdb/qbb2020-answers/day4-lunch/", "BDGP6.Ensembl.81.gtf"))

get_the_output(21378950, read_gtf_df("/Users/cmdb/qbb2020-answers/day4-lunch/", "BDGP6.Ensembl.81.gtf"), closest_index)

#Output:
#------------------------------------------------------------------------------------------------------
#This thing iterated 17 times before finding the answer
#The gene name is tin
#The gene id is FBgn0004110
#The biotype is protein_coding
#Shortest distance from the beginning of gene is 27 bp
#------------------------------------------------------------------------------------------------------

