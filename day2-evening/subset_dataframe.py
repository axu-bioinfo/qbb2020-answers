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

common_variation_final = (alt_allele_freq_SNP_df > 0.05) & (alt_allele_freq_SNP_df < 0.95)
common_variation_sub_df = alt_allele_freq_SNP_df.loc[common_variation_final]
print(common_variation_sub_df)

#did it again but I think the logic is pretty cool
common_variation_part1 = alt_allele_freq_SNP_df > 0.05
common_variation_part2 = alt_allele_freq_SNP_df < 0.95

common_variation_part3 = pd.concat([common_variation_part1, common_variation_part2], axis = 1)
common_variation = common_variation_part3.prod(axis = 1)

common_variation = common_variation.astype(bool)
common_variation_sub_df = alt_allele_freq_SNP_df.loc[common_variation]
print(common_variation_sub_df)

