import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#count where to start reading vcf data
count = 0
with open ("/Users/cmdb/Desktop/quant_assignment_2/freebayes_output/annotated_func_A01_sacCer3_variants.vcf", "r") as refFile:
	lines = refFile.readlines()
	for line in lines:
		if "##" in line:
			count = count + 1
		if "##" not in line:
			break

#obtain the format and info columns
vcf_df = pd.read_csv("/Users/cmdb/Desktop/quant_assignment_2/freebayes_output/annotated_func_A01_sacCer3_variants.vcf", skiprows=count + 1, sep="\t", header=None)
format_series = vcf_df[7]
info_series = vcf_df[7]

depth_distribution = []
quality_distribution = []
allele_freq = []
loss_percentage_effected = []
for index, value in format_series.items():
	DP_start = value.find("DP=") + 3
	DP_end = value.find(";", DP_start)
	depth_distribution.append(value[DP_start:DP_end])

	QA_start = value.find("QA=") + 3
	QA_end = value.find(";", QA_start)
	quality_distribution.append(value[QA_start:QA_end])

for index, value in info_series.items():
	AF_start = value.find(";AF=") + 4
	AF_end = value.find(";", AF_start)
	allele_freq.append(value[AF_start:AF_end])

	if "LOF=" in value:
		LOF_start = value.find("LOF=") + 4
		LOF_end = value.find(")", LOF_start)
		LOF_segment = value[LOF_start:LOF_end]
		LOF_value = LOF_segment[LOF_segment.rfind("|") + 1:]
		loss_percentage_effected.append(LOF_value)

fig, ax = plt.subplots(2, 2)
ax0, ax1, ax2, ax3 = ax.flatten()

ax0.hist(depth_distribution)
ax0.title.set_text("Read depth distribution")
ax0.set_xlabel("Read depth of variant genotype")
ax0.set_ylabel("Number of times the read depth has been reached for total sampled variants")

ax1.hist(quality_distribution)
ax1.title.set_text("Quality distribution")
ax1.set_xlabel("Quality of variant genotype")
ax1.set_ylabel("Number of times the quality has been reached for total sampled variants")

ax2.hist(allele_freq)
ax2.title.set_text("Allele frequency spectrum of variants")
ax2.set_xlabel("Allele frequency of variant")
ax2.set_ylabel("Number of times the allele frequency has been reached for total sampled variants")

ax3.hist(loss_percentage_effected)
ax3.title.set_text("Summary of loss of function prediction percent affected by variant")
ax3.set_xlabel("Loss of function predicted percent affected")
ax3.set_ylabel("Number of times the percent affected has been reached for total sampled variants")

plt.savefig("week2_subplots.png")
plt.show()
