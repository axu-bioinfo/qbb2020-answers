# PCA of SNP data

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

df_matrix_1kg = pd.read_csv("/Users/cmdb/Desktop/matrix_1kg.txt" , sep ='\t')

#obtain list of samples and use it to get a sub df of the original data
list_of_sample = df_matrix_1kg.columns.values.tolist()
all_SNP_count_df = df_matrix_1kg[list_of_sample[4:]]

#add the SNP counts together
SNP_sum_df = all_SNP_count_df.sum(axis = 1)
total_samples = len(list_of_sample[4:])
total_chromosomes = total_samples * 2

#divide number of chromosomes by total count of SNP variations within population
alt_allele_freq_SNP_df = SNP_sum_df.div(total_chromosomes)

#plot the histogram
fig, ax = plt.subplots()
ax.bar(alt_allele_freq_SNP_df.index, alt_allele_freq_SNP_df)
ax.title.set_text("SNP alt frequency")
ax.set_xlim([0.0, 10005])
ax.set_ylim([0, 1])
ax.set_xlabel("SNP index")
ax.set_ylabel("frequency")
plt.show()