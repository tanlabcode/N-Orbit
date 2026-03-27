import pandas as pd
import glob
import numpy as np

''' This script compiles the neighborhood distance columns generated in Step 1a into the full neigborhood distace matrix.'''
# Input paths and column name settings (same as Step 1)
input_file_path = "data/synthetic_mrf_neighborhoods_v1.csv"
intermediate_path = "intermediates/SyntheticV1/"
im_label = "Image"
neighborhood_label = "Neighborhood"

cells = pd.read_csv(input_file_path)
norbits = pd.read_csv(intermediate_path+ "norbits.csv")

# Compile neighborhood distance vectors into matrix
neighborhood_list = sorted(set(norbits["unit"]))
dfs = []
for file in glob.glob(intermediate_path + "dists/*.csv"):
    dfs.append(pd.read_csv(file).drop(["Unnamed: 0"], axis=1))
df = pd.concat(dfs, axis = 1)
df = df.loc[:,~df.columns.duplicated()].copy()
df = df.reindex(sorted(df.columns), axis=1)
df = pd.DataFrame(df.values+np.transpose(df.values), index=list(df.columns),columns=list(df.columns))
df.to_csv(intermediate_path+"neighborhood_distance_matrix.csv")
