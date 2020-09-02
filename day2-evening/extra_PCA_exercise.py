# The extra exercise

import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def obtain_PCA_of_SNP(reference_text_file):
	df_matrix_1kg = pd.read_csv(reference_text_file , sep ='\t')

	#obtain list of samples and use it to get a sub df of the original data
	list_of_sample = df_matrix_1kg.columns.values.tolist()
	all_SNP_count_df = df_matrix_1kg[list_of_sample[4:]]

	#add the SNP counts together
	SNP_sum_df = all_SNP_count_df.sum(axis = 1)
	total_samples = len(list_of_sample[4:])
	total_chromosomes = total_samples * 2

	#divide number of chromosomes by total count of SNP variations within population
	alt_allele_freq_SNP_df = SNP_sum_df.div(total_chromosomes)

	#get the common variation "mask" to find the sample dataframes (not added together)
	common_variation_final = (alt_allele_freq_SNP_df > 0.05) & (alt_allele_freq_SNP_df < 0.95)
	common_variation_pca_df = all_SNP_count_df.loc[common_variation_final]

	#samples as rows and the SNP alt as column
	common_variation_pca_T_df = common_variation_pca_df.T

	#standardize the data
	standardized_SNP = StandardScaler().fit_transform(common_variation_pca_T_df)
	#print(standardized_SNP.mean(axis = 0))
	#print(standardized_SNP.var(axis = 0))

	#find the top 5 principle components
	pca = PCA(n_components = 5)
	pca_output = pca.fit_transform(standardized_SNP)
	pca_output

	#convert the pca output to pandas
	pca_output_df = pd.DataFrame(data = pca_output, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
	pca_output_df['sample'] = common_variation_pca_T_df.index
	return(pca_output_df)

def add_metadata(metadata_file, parameter):
	metadata_df = pd.read_csv(metadata_file , sep ='\t')
	metadata_df = metadata_df[["sample", parameter]]
	return(metadata_df)


#population
PCA_SNP_w_label = obtain_PCA_of_SNP("/Users/cmdb/Desktop/matrix_1kg.txt")
pop_df = add_metadata("/Users/cmdb/Desktop/integrated_call_samples_v3.20130502.ALL.panel", "pop")

pop_PCA = pd.merge(PCA_SNP_w_label, pop_df, on = "sample")

#plot the dataframe with tissue color as labeled in the metadata
fig, ax = plt.subplots()
#add color
population = pop_PCA.groupby("pop")
for name, group in population:
    ax.scatter(x = group['PC1'], y = group['PC2'], label = name)
ax.title.set_text("PCA of 2548 people's SNPs colored by country")
plt.legend(bbox_to_anchor=(1.04, 1), loc = 'lower left', ncol = 4)
plt.show()

#superpopulation

#sex




