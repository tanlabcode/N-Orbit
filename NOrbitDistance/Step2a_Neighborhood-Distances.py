import pandas as pd
import numpy as np
from iteround import saferound
from scipy.spatial import distance_matrix, minkowski_distance
from munkres import Munkres
import plotly.express as px
from matplotlib import pyplot as plt
import sys
import os
from scipy.optimize import linear_sum_assignment

''' This script calculates neighborhood distances for all neighborhoods specified in UNIT to all neighborhoods in the dataset that precede it alphabetically.
'''

# If UNIT is image, use UNIT_MODE = "Image". If UNIT is neighborhood, use UNIT_MODE = "Neighborhood"
# To run all images in a single job, use UNIT_MODE = "All" (slow)
UNIT = sys.argv[1]
UNIT_MODE = "All"

# Path to intermediates (same as Step1)
intermediate_path = "intermediates/SyntheticV1/"

# Read N-Orbits from Step1 and gather list of all neighborhoods/units in dataset
df = pd.read_csv(intermediate_path + "norbits.csv")
neighborhood_list = sorted(set(df["unit"]))

# Filter neighborhood list to that specified in command line for parallel processing.
if UNIT_MODE == "Image":
	unit_neighborhoods = [item for item in neighborhood_list if item.startswith(UNIT)]
elif UNIT_MODE == "Neighborhood":
	unit_neighborhoods = [item for item in neighborhood_list if item == UNIT]
elif UNIT_MODE == "All":
	unit_neighborhoods = list(neighborhood_list)

# Create directory if not already exists
if not os.path.exists(intermediate_path+"dists/"):
    os.makedirs(intermediate_path+"dists/")

print(unit_neighborhoods)

def neighborhood_distance(df1_vectors, df2_vectors):
    '''Get neighborhood distance between a pair of neighborhoods.
    df1_vectors: Pandas DataFrame of sampled N-Orbit vectors for Neighborhood 1
    df2_vectors: Pandas DataFrame of sampled N-Orbit vectors for Neighborhood 2
    Output: Pairwise distance between Neighborhood 1 and 2, positions for linear sum assignment (not used)
    '''
    distance_mat = distance_matrix(df1_vectors, df2_vectors, p = 1)
    distance_mat_orig = np.copy(distance_mat)
    row_ind, column_ind = linear_sum_assignment(distance_mat)
    return np.sum(distance_mat[row_ind,column_ind]), row_ind, column_ind

import warnings
warnings.filterwarnings("ignore")

def main():
    unit_dists = {}
    
    n_vectors = {}
    print("Starting vector writing.")
    # Get NumPy formulation of N-Orbit vectors for all neighborhoods.
    for i in neighborhood_list:
            df1_vectors = pd.read_csv(intermediate_path + "sampled_vectors/"+i+".csv").set_index("Unnamed: 0")
            n_vectors[i] = df1_vectors.values
    print("Finished vector writing.")
    
    # For all specified neighborhoods in this job, calculate neighborhood distances to all neighborhoods in the dataset that come before it alphabetically.
    for i in unit_neighborhoods:
        unit_dists[i] = np.zeros((len(neighborhood_list)))
        print(i)
        for j in range(len(neighborhood_list)):
            if i > neighborhood_list[j]:
                dist, row_ind, column_ind = neighborhood_distance(n_vectors[i],n_vectors[neighborhood_list[j]]) 
                unit_dists[i][j] = dist
            else:
                break
                
    # Save columns of distance matrix for neighborhoods specified for this job
    unit_dists = pd.DataFrame(unit_dists)
    unit_dists.to_csv(intermediate_path+"dists/" + UNIT+"_ndists.csv")

if __name__ == "__main__":
    main()
