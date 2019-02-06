import random
import pandas as pd

pca_dataframe = pd.DataFrame()

for pc in ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']:

    sample_headers = []
    sample_values = []

    for sample in range(1, 1001):

        sample_headers.append('Sample_%s' % sample )
        sample_values.append(round(random.uniform(-10, 10),4))

    pc_series = pd.Series(sample_values, index = sample_headers)


    pca_dataframe[pc] = pc_series

pca_dataframe.to_csv('PCA_Data.tsv', sep = '\t')
