#!/usr/bin/env python3

import os, sys
from fasta_iterator_class import FASTAReader


def obtain_shell_input(input):
	return(input)

shell_input = obtain_shell_input(sys.argv)

target = sys.argv[1]
query = sys.argv[2]
k = sys.argv[3]

def get_kmers(reader, k):
	kmers = {}
	k = int(k)
	for seq_id, sequence in reader:
	    for i in range(0, len(sequence) - k+1):
	        kmer = sequence[i:i + k]
	        kmers.setdefault(kmer, 0)
	        kmers[kmer] =[seq_id, i]
	return(kmers)

target_reader = FASTAReader(open(target))
target_dic = get_kmers(target_reader, k)

query_reader = FASTAReader(open(query))
query_dic = get_kmers(query_reader, k)

print("target sequence name\ttarget_start\tquery_start\tkmer")

for query_key in query_dic:
	for target_key in target_dic:
		if query_key == target_key:
			print(target_dic[query_key][0]  + "\t" + str(target_dic[query_key][1]) + "\t" + str(query_dic[query_key][1]) + "\t" + query_key)

if __name__ == "__main__":
	True