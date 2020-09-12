#!/usr/bin/env python3

import os, sys
from fasta_iterator_class import FASTAReader

#run this on shell
# >> ./kmer_matcher.py subset.fa droYak2_seq.fa 11

def obtain_shell_input(input):
	return(input)

shell_input = obtain_shell_input(sys.argv)

target = sys.argv[1]
query = sys.argv[2]
k = sys.argv[3]

#function to get a dictionary of kmers
#kmer sequence as key, and a list of lists containing the sequence id and the kmer location
#on the fasta sequence the kmer came from
#if there is more than one instance of a kmer it will append another list of [seq_id, location]
#to the dictionary kmer sequence list value
def get_kmers(reader, k):
	kmers = {}
	k = int(k)
	for seq_id, sequence in reader:
	    for i in range(0, len(sequence) - k+1):
	    	kmer = sequence[i:i + k]
	    	if kmer in kmers:
	    		kmers[kmer].append([seq_id, i])
	    	else:
		        kmers[kmer] =[[seq_id, i]]
	return(kmers)

target_reader = FASTAReader(open(target))
target_dic = get_kmers(target_reader, k)

query_reader = FASTAReader(open(query))
query_dic = get_kmers(query_reader, k)

with open("kmer_matcher_first1000.out", "a+") as newFile:
	newFile.write("target_sequence_name\ttarget_start\tquery_start\tk-mer\n")

count = 0
for query_key in query_dic:
	for target_key in target_dic:
		if query_key == target_key:
			if count > 999:
				sys.exit()
			if len(target_dic[query_key]) == 1:
				count = count + 1
				with open("kmer_matcher_first1000.out", "a+") as newFile:
					newFile.write(target_dic[query_key][0][0]  + "\t" + str(target_dic[query_key][0][1]) + "\t" + str(query_dic[query_key][0][1]) + "\t" + query_key + "\n")
			else:
				for num_kmer_repeat in range(len(target_dic[query_key])):
					count = count + 1
					with open("kmer_matcher_first1000.out", "a+") as newFile:
						newFile.write(target_dic[query_key][num_kmer_repeat][0]  + "\t" + str(target_dic[query_key][num_kmer_repeat][1]) + "\t" + str(query_dic[query_key][0][1]) + "\t" + query_key + "\n")

if __name__ == "__main__":
	True