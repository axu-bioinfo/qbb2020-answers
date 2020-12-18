import os, sys
import pandas as pd
import re

def bin_kraken(bin_file, kraken_df):
	kraken_node_list = kraken_df[0].to_list()
	bin_df = pd.read_csv(bin_file, sep = ">", header = None)
	bin_node_list = bin_df[1].to_list()
	for node in bin_node_list:
		if node not in kraken_node_list:
			continue
		index = kraken_node_list.index(node)
		taxon = kraken_df.iloc[index][1]
		taxon_tab = re.sub(";", "\t", taxon)
		with open(bin_file + ".kraken", "a+") as newFile:
			newFile.write(taxon_tab + "\n")


kraken_df = pd.read_csv("assembly.kraken", sep = "\t", header = None)

#bin_kraken("bin.1.fa.node.name", kraken_df)
#bin_kraken("bin.2.fa.node.name", kraken_df)
#bin_kraken("bin.3.fa.node.name", kraken_df)
#bin_kraken("bin.4.fa.node.name", kraken_df)
#bin_kraken("bin.5.fa.node.name", kraken_df)
#bin_kraken("bin.6.fa.node.name", kraken_df)
#bin_kraken("bin.7.fa.node.name", kraken_df)
#bin_kraken("bin.8.fa.node.name", kraken_df)

for i in range(8):
	print("ktImportText -q bin." + str(i+1) + ".fa.node.name.kraken -o /Users/cmdb/Desktop/week10_test_html/bin_html/bin." + str(i+1) + ".html")
