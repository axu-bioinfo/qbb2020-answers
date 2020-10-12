import os, sys
from Bio import SeqIO

#convert fasta file to dictionary
#the dictionary has the sequence id as the key and the sequence as the value
def fasta_file_to_dic(file):
	fasta_dic = {}
	for seq_record in SeqIO.parse(file, "fasta"):
		if seq_record.id in fasta_dic:
			print("Uh oh!")
		fasta_dic[seq_record.id] = seq_record.seq
	return(fasta_dic)

nucleotide_dic = fasta_file_to_dic("/Users/cmdb/Desktop/quant_hw4/final_seq.fasta")
protein_dic = fasta_file_to_dic("/Users/cmdb/Desktop/quant_hw4/out_seq.fasta")

#match aa to nucleotide which prints a dictionary prot_nuc_match_dic
#the key value is the protein fasta file id and the value is a list of lists
#each value has the amino acid and the corresponding nucleotide
def prot_nuc_match_dic(nuc_dic, prot_dic):
	prot_nuc_match_dic = {}
	for protein_id in prot_dic:
		nucleotide_id = protein_id[ : protein_id.rfind("_")]
		protein = prot_dic[protein_id]
		nucleotide = nuc_dic[nucleotide_id]
		nuc_length = 0
		aa_and_nuc_list = []
		for aa in protein:
			codon = nucleotide[nuc_length : nuc_length + 3]
			nuc_length = nuc_length + 3
			aa_and_nuc_list.append([aa, str(codon)])
		prot_nuc_match_dic[protein_id] = aa_and_nuc_list
	return(prot_nuc_match_dic)

prot_nuc_dic = prot_nuc_match_dic(nucleotide_dic, protein_dic)

#go through the mafft file and get a dic of dic[prot_id] = [fasta description, alignment with gaps]
def mafft_to_dic(file):
	fasta_dic = {}
	for seq_record in SeqIO.parse(file, "fasta"):
		if seq_record.id in fasta_dic:
			print("Uh oh!")
		fasta_dic[seq_record.id] = [seq_record.description, str(seq_record.seq)]
	return(fasta_dic)

mafft_dic = mafft_to_dic("/Users/cmdb/Desktop/quant_hw4/mafft_align_out.fasta")

def mafft_nuc(mafft_dic, prot_nuc_dic):
	mafft_nuc_dic = {}
	for prot_id in mafft_dic:
		mafft_nuc = ""
		prot_seq_index = 0
		for alignment in mafft_dic[prot_id][1]:
			if alignment == "-":
				mafft_nuc = mafft_nuc + "---"
			else:
				mafft_nuc = mafft_nuc + prot_nuc_dic[prot_id][prot_seq_index][1]
				if prot_seq_index + 1 < len(prot_nuc_dic[prot_id]) and prot_nuc_dic[prot_id][prot_seq_index + 1][0] == "*":
					prot_seq_index = prot_seq_index + 1
				prot_seq_index = prot_seq_index + 1
		mafft_nuc_dic[mafft_dic[prot_id][0]] = mafft_nuc
	return(mafft_nuc_dic)

mafft_nuc_dic = mafft_nuc(mafft_dic, prot_nuc_dic)

for header in mafft_nuc_dic:
	with open("mafft_nuc_aligned.txt", "a+") as newFile:
		newFile.write(">" + header + "\n" + mafft_nuc_dic[header] + "\n")

os.system("mv mafft_nuc_aligned.txt mafft_nuc_aligned.fasta")


