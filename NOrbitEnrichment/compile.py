import pandas as pd
import glob
import numpy as np
from collections import Counter
from scipy import stats

# Set paths and column labels
trial_path = "../trials/"
input_file_path = "../data/synthetic_mrf_neighborhoods_v1.csv"
cell_type_label = "Cell_Type"
output_path = "../to/output/"

# Compile trials into single dataframe
dfs = []
trials = glob.glob(trial_path+"*.csv")
for f in trials:
    dfs.append(pd.read_csv(f))
df = pd.concat(dfs)
df = df.drop(columns = "Unnamed: 0")

# Recalculate p and q values
cellTypes = sorted(set(pd.read_csv(input_file_path)[cell_type_label]))
df = pd.DataFrame(df.groupby(list(df.columns[:2*len(cellTypes)])).mean()["pvalue"]).reset_index().sort_values("pvalue")
df["qvalue"] = stats.false_discovery_control(df["pvalue"])
df = df.reset_index(drop=True)
df = df.sort_values("pvalue")
df = df.reset_index(drop=True)
df.to_csv(output_path+"n-orbit-enrichment-results.csv")
