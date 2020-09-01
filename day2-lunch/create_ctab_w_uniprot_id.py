# rewrite the t_data.ctab file so that the FlyBase gene ID is replaced with the UniProt ID

# Directions: in the last line, change the text file locations to the input file, mapping file, and output file
# and then run the python script
import pandas as pd

def rewrite_the_ctab(input_file, mapping_file, output_file):
	mapping_file_dic = {}
	mapping_file_dic = dic_mapping_file(mapping_file)
	with open(input_file, "r") as refFile:
		all_lines = refFile.readlines()
		count = 0
		with open(output_file, "a+") as final_outFile:
			final_outFile.write(all_lines[0])
		for line in all_lines[1:]:
			new_string = ""
			if count > 99:
				break
			line = line.split('\t')
			flybase_gene_id_ctab = line[8]
			if flybase_gene_id_ctab in mapping_file_dic:
				line[8] = mapping_file_dic[flybase_gene_id_ctab]
			else:
				line[8] = "*"
			for char in line:
				if "\n" in char:
					new_string = new_string + char
				else:
					new_string = new_string + char + "\t"
			with open(output_file, "a+") as final_outFile:
				final_outFile.write(new_string)
			count = count + 1


def dic_mapping_file(mapping_file):
	mapping_file_dic = {}
	mapping_file_df = pd.read_csv(mapping_file, sep='\t', header=None)
	for val in mapping_file_df.values:
		if val[0] == "NO GENE NAME":
			continue
		mapping_file_dic[val[0]] = val[1]
	return(mapping_file_dic)


rewrite_the_ctab("/Users/cmdb/Desktop/t_data.ctab", "/Users/cmdb/Desktop/mapping_file.txt", "/Users/cmdb/Desktop/uniprot_t_data.ctab")