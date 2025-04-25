# In-house Python codes used in Chapter 3: core genome LIN code for H. influenzae

import pandas as pd
import os
import numpy as np

## 1) convert distance matrix of pairwise allelic mismatch to frequency table

os.chdir('dir/to/input_file')

infile = 'distance_matrix.xlsx' # distance matrix from PubMLST Genome Comparator
df_matrix = pd.read_excel(infile, sheet_name = 'distance matrix', index_col=0)

df_freq = df_matrix.stack().rename_axis(['genome1','genome2']).reset_index(name='allelic_mismatch')
df_freq['genome_combine'] = df_freq['genome1'].astype(str) + df_freq['genome2'].astype(str)
df_freq.set_index('genome_combine', inplace = True)

## 2) k-means clustering

from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

lst_allelic_mismatch = df_freq['allelic_mismatch'].to_list()
arr_allelic_mismatch = np.array(lst_allelic_mismatch)
arr_allelic_mismatch = arr_allelic_mismatch.reshape(-1, 1) #reshape for one-dimensional data

inertias = []
for i in range(2,8): #lowest number cannot be 1, highest number was chosen arbitrarily
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(arr_allelic_mismatch)
    inertias.append(kmeans.inertia_)

plt.plot(range(2,8), inertias, marker = 'o')
plt.show() # visualise to choose N clusters based on drop in intertia

kmeans = KMeans(n_clusters=5, init='k-means++', random_state=100)
kmeans.fit(arr_allelic_mismatch)

clusters = []
clusters = kmeans.labels_

df_freq['cluster'] = clusters

outfile = 'frequency_table_kmeans_cluster.csv'

df_freq.to_csv(outfile)

### --- ###


## 3) Rand index

from sklearn.metrics.cluster import adjusted_rand_score

os.chdir('dir/to/input_file')

infile = 'dataset_all_clustering_results.csv' # summary of clustering result for each genome in the dataset using different clustering methods
df_clust = pd.read_csv(infile)

non_lin = list(range(0,9))
lin = list(range(9,17))

list1 = []
list2 = []
list3 = []

for num in non_lin:
    for no in lin:
        df_tmp = df_clust.iloc[: , [num,no]].copy()
        df_tmp = df_tmp.dropna()
        headers = list(df_tmp)
        list1.append(headers[0])
        list2.append(headers[1])
        non_lin_list = df_tmp[df_tmp.columns[0]].values.tolist()
        lin_list = df_tmp[df_tmp.columns[1]].values.tolist()
        #lincode cluster as the "ground truth class labels"
        score = adjusted_rand_score(lin_list, non_lin_list)
        list3.append(score)

dfresult2 = pd.DataFrame(
        {
            'Non-LINcode clustering' : list1,
            'LINcode clustering threshold' : list2,
            'Adjusted rand idx' : list3
        }
    )

outfile = 'rand_index_cluster_comparison.csv'
dfresult2.to_csv(outfile)
