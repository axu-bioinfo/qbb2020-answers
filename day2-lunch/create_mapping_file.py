# python script to make a mapping file of FlyBase IDs and UniProt IDs

# Directions: run this script by changing the location of the fly.txt text file and where you want
# the output to be created, and run the python script

def create_mapping_file(where_fly_txt, where_output_file):
	with open(where_fly_txt, "r") as refFile:
		eachlines = refFile.readlines()
		for line in eachlines:
			parsed_line = []
			if "DROME" in line:
				line = line.split(" ")
				for item in line:
					if item != '':
						parsed_line.append(item)
				#the genes with no name are renamed "NULL"
				with open(where_output_file, "a+") as newFile:
					if "\n" == parsed_line[-1]:
						parsed_line[-1] = "NO GENE NAME"
					newFile.write(parsed_line[-1].strip() + "\t" + parsed_line[-2] + "\n")


# should output a tab delimited text file contianing first the FlyBase ID and then the UniProt ID
create_mapping_file("/Users/cmdb/Desktop/fly.txt", "/Users/cmdb/Desktop/mapping_file.txt")