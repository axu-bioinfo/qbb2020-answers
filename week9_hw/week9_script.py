#!/usr/bin/env python2

import numpy as np
import hifive
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pyBigWig

def heatmap():
	hic = hifive.HiC("genome/filter_normal.h5", "r")

	data = hic.cis_heatmap('chr13', 1000000, datatype = "fend", arraytype = "full", diagonalincluded = True)

	enrichment = (data[:, :, 0] + 1) / (data[:, :, 1] + 1)

	enrichment = np.log(enrichment)

	fig, ax = plt.subplots()
	c = ax.pcolor(enrichment)
	ax.set_title("Heatmap for Chromosome 13")
	plt.show()

def make_bed():
	hic = hifive.HiC("genome/filter_normal.h5", "r")
	Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
	Comp.write_eigen_scores('hic_comp.bed')

def get_comp_scores():
	hic = hifive.HiC("genome/filter_normal.h5", "r")
	Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
	X = Comp.positions['chr13']
	Y = Comp.eigenv['chr13']
	new_X = [gene_loc[0] for gene_loc in X]
	fig, ax = plt.subplots()
	ax.step(new_X, Y)
	ax.set_title("Compartment scores")
	ax.set_xlabel("Compartment Score")
	ax.set_ylabel("Genomic Position")
	plt.show()

# to get the positive and negative values separated
# >> grep "-" hic_comp.bed > neg_hic_comp.bed
# >> grep -i -v -E "-" hic_comp.bed > pos_hic_comp.bed
# >> bedtools intersect -a data/WT_fpkm.bed -b neg_hic_comp.bed > neg_intersect.bed
# >> bedtools intersect -a data/WT_fpkm.bed -b pos_hic_comp.bed > pos_intersect.bed

def get_last_col(filename):
	bed = pd.read_csv(filename, sep = "\t", header = None)
	bed_list = bed[4].to_list()
	return(bed_list)

def get_col(filename, col):
	bed = pd.read_csv(filename, sep = "\t", header = None)
	bed_list = bed[col].to_list()
	return(bed_list)

def gene_expression_violin_plot():
	pos_val = get_last_col("pos_intersect.bed")
	neg_val = get_last_col("neg_intersect.bed")
	data = [neg_val, pos_val]
	fig, ax = plt.subplots()
	ax.set_ylabel("FPKM")
	ax.set_xlabel("Negative (left) Positive (right)")
	ax.violinplot(data)
	ax.set_title("Violin plot for gene expression")
	plt.show()

def neg_methylation_expression():
	neg_start = get_col("neg_intersect.bed", 1)
	neg_end = get_col("neg_intersect.bed", 2)
	bw = pyBigWig.open('data/WT_H3K27me3.bw')
	neg_methyl = []
	for i in range(len(neg_start)):
		stat = bw.stats('chr13', neg_start[i], neg_end[i], type='sum')
		if stat[0] == "None":
			stat[0] = 0
		neg_methyl.append(stat[0])
	neg_expression = get_last_col("neg_intersect.bed")

	fig, ax = plt.subplots()
	ax.set_ylabel("H3K27me3")
	ax.set_xlabel("Expression")
	ax.scatter(neg_expression, neg_methyl, alpha = 0.5)
	ax.set_title("A compartment H3K27me3 vs Expression")
	plt.show()

def pos_methylation_expression():
	pos_start = get_col("pos_intersect.bed", 1)
	pos_end = get_col("pos_intersect.bed", 2)
	bw = pyBigWig.open('data/WT_H3K27me3.bw')
	pos_methyl = []
	for i in range(len(pos_start)):
		stat = bw.stats('chr13', pos_start[i], pos_end[i], type='sum')
		if stat[0] == "None":
			stat[0] = 0
		pos_methyl.append(stat[0])
	pos_expression = get_last_col("pos_intersect.bed")

	fig, ax = plt.subplots()
	ax.set_ylabel("H3K27me3")
	ax.set_xlabel("Expression")
	ax.scatter(pos_expression, pos_methyl, alpha = 0.5)
	ax.set_title("B compartment H3K27me3 vs Expression")
	plt.show()


