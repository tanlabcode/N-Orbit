import numpy as np
import pandas as pd
from permutation_test import enrichment_analysis
from scipy import stats
import sys
import os

vectors_path = "/path/to/cluster/vectors/NSCLC-D-cluster-0.csv"
output_path = "/path/to/output/"
run = sys.argv[1]

if os.path.isdir(output_path):
    print(f"'{output_path}' is a directory.")
else:
    os.mkdir(output_path)

# Number of permutations
# If doing multiple runs, num_permutations needs to be the same for all runs
num_permutations = 25000

# Number of cell types (first numCellTypes columns are not permuted)
numCellTypes = 16

# Threshold for binarization (equivalent to 50um)
threshold = 0.2

# Load data
tcn_vectors_df = pd.read_csv(vectors_path)
tcn_vectors_df = tcn_vectors_df.drop(columns="neighborhood")
tcn_vectors = (tcn_vectors_df.values>threshold).astype(np.int32)
tcn_vectors_df = pd.DataFrame(tcn_vectors)

# Subsample size
sample_size = int(round(len(tcn_vectors_df)*0.2))

# Perform the enrichment analysis
p_values = enrichment_analysis(tcn_vectors, num_permutations, numCellTypes, sample_size)

# Add p-values to the original dataframe
tcn_vectors_df['pvalue'] = p_values
tcn_vectors_df = tcn_vectors_df.drop_duplicates()
tcn_vectors_df = tcn_vectors_df.sort_values("pvalue")
tcn_vectors_df['qvalue'] = stats.false_discovery_control(tcn_vectors_df['pvalue'], method='bh')

# Save the results to a CSV file
tcn_vectors_df.reset_index(drop=True).to_csv(output_path + "NOrbit_Enrichment_Run_"+str(run)+".csv")
