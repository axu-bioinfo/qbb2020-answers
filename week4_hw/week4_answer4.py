import os, sys
import pandas as pd
import numpy as np
from Bio import SeqIO

#count the dN and dS, if the dna sequence doesn't match but the aa does then it is a dS
#if the dna and the aa doesn't match, it counts as a dN

#first obtain a dataframe of all the aligned nucleotide codons and proteins
def fasta_to_df(fasta_file, n):
	out_dic = {}
	for seq_record in SeqIO.parse(fasta_file, "fasta"):
		out_list = []
		out_count = 0
		while out_count < len(str(seq_record.seq)):
			out = str(seq_record.seq)[out_count : out_count + n]
			out_count = out_count + n
			out_list.append(out)
		out_dic[seq_record.id] = out_list

	out_df = pd.DataFrame(out_dic)
	return(out_df)

#codon table from the assignment
codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

#calculate the dN and dS of each column
def get_dN_dS(codon_df, protein_df):
	dn_ds_dic = {}
	codon_trans_df = codon_df.T
	protein_trans_df = protein_df.T
	for i in codon_trans_df.columns:
		#1 means codon matches, 0 means codon doesn't match, first value will always be 1
		codon_bool = np.where(codon_trans_df[i]["query_1"] == codon_trans_df[i], 1, 0)

		#1 means codon matches, 0 means codon doesn't match
		protein_bool = np.where(protein_trans_df[i]["query_1"] == protein_trans_df[i], 1, 0)

		#add them together, if it equals 2 then both codon + aa were the same
		#if it equals 1, then the dna sequence doesn't match, but the aa does --> dS
		#if it equals 0, then both the protein and the dna sequence don't match --> dN
		dN_dS_bool = codon_bool + protein_bool
		dS = np.count_nonzero(dN_dS_bool == 1)
		dN = np.count_nonzero(dN_dS_bool == 0)
		d_list = [dS, dN]
		#i stands for the codon position
		dn_ds_dic[i] = d_list

	#for the return pd, the 0 row is the dS value, 1 row is the dN value
	dn_ds_df = pd.DataFrame(dn_ds_dic)
	return(dn_ds_df)

#use the aligned codon file name
codon_df = fasta_to_df("mafft_nuc_aligned.fasta", 3)
#change the file name to the protein aligned file name
protein_df = fasta_to_df("mafft_align_out.fasta", 1)

dN_dS_df = get_dN_dS(codon_df, protein_df)
dN_dS_df.to_csv("ds_dn.txt", sep = "\t", index = None)




