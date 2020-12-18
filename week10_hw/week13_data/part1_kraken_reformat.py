import os, sys
import pandas as pd
import re

def convert_file(filename):
	#read the original kraken file to a dataframe
	original_df = pd.read_csv(filename, sep = "\t", header = None)
	#check for duplicate samples and put the sample id in a dic with sampled id as key and
	#[sample count, heirarchies] as values
	id_name_list = original_df[0].to_list()
	heirarchies_dic = {}
	for i in range(len(id_name_list)):
		if original_df.iloc[i][0] in heirarchies_dic:
			heirarchies_dic[original_df.iloc[i][0]][0] = heirarchies_dic[original_df.iloc[i][0]][0] + 1
		else:
			heirarchies_dic[original_df.iloc[i][0]] = [1, original_df.iloc[i][1]]
	
	for key in heirarchies_dic:
		count = heirarchies_dic[key][0]
		heirarchy = re.sub(";", "\t", heirarchies_dic[key][1])
		with open(filename + ".mod", "a+") as newFile:
			newFile.write(str(count) + "\t" + heirarchy + "\n")
	

#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492183.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492186.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492188.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492189.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492190.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492193.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492194.kraken")
#convert_file("/Users/cmdb/qbb2020-answers/week10_hw/week13_data/KRAKEN/SRR492197.kraken")