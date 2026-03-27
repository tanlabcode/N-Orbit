![header](https://github.com/barbara-xiong/N-Orbit/blob/main/images/NOrbitLogo-03.png)

# **Distance-Based Neighborhood Exploration with N-Orbits**

## **Contents**

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Maintainers](#maintainers)
- [Citation](#citation)

## **Overview**

![N-Orbit Schematic](https://github.com/barbara-xiong/N-Orbit/blob/main/images/Fig1_NOrbitSchematic.png)

Tissue cellular neighborhoods (TCNs) are spatially contiguous regions of homogeneous and distinct cell type composition. Formation of TCNs may be indicative of various cell types coordinating to form functional niches. Studying spatial relationships on the level of TCNs may allow for identification of localized patterns that may otherwise be diluted and overlooked when analyzing in bulk.

Several neighborhood detection methods have been developed in recent years including our lab's own method, [CytoCommunity](https://github.com/huBioinfo/CytoCommunity). However, downstream methods to quantify TCN changes across conditions (e.g. time, clinical subtypes) beyond the level of cell type enrichment are limited. To achieve this goal, we introduce a distance-based approach, centered around a novel N-Orbit formalism of a neighborhood.

The N-Orbit structure, encodes the proximity of each possible cell type to each cell in a TCN. The focus on proximity, rather than exact subgraph connections, allows for increased computational efficiency and greater instances of each structure. N-Orbits can be formulated as vectors, allowing for calculation of pairwise distance between two N-Orbit structures using their Manhattan distance. Distances between neighborhoods can then be computed from the minimum total distance between their representative set of N-Orbits. Pairwise neighborhood distances can be compiled into an overall neighborhood distance matrix that can be used for clustering, projection, visualization, and other methods for analyzing neighborhoods on a global scale.

## **Installation**

### **Hardware requirement**

CPU: i7

Memory: 16G or more

Storage: 10G or more

### **Software requirement**

Conda version: 4.13.0

Python version: 3.10.4

Clone this repository and cd into it as below

```bash
git clone https://github.com/barbara-xiong/N-Orbit.git
cd N-Orbit
```

**For Linux**

**Preparing the virtual environment**

```bash
conda env create -f environment_norbit.yml
conda activate NOrbit
```

**For Windows**

**Preparing the virtual environment**

```bash
conda env create -f environment_norbit_win.yml
conda activate NOrbit
```

## **Usage**

This package requires single-cell resolution spatial omics data, where cell types have already been annotated.
By default, the N-Orbit distance pipeline is intended for neighborhood-level use. However, small adjustments may be made to calculate distances on the sample level and are italicized below.

**Prepare input data**

To prepare the input data, compile a CSV table where each row represents a cell and columns are provided for **(1) x coordinates, (2) y coordinates, (3) cell type labels, (4) sample/image labels, and (5) neighborhood/TCN labels**. Column names are flexible and specified in Step1 below. Synthetic data examples are provided under the directory examples/synthetic_mrf_neighborhoods_v1.csv and examples/synthetic_mrf_neighborhoods_v2.csv.

*For calculating sample-level distances, use a column of a constant values (e.g. all zeros) for the neighborhood label column.*

### **N-Orbit Distance**

Scripts for calculating N-Orbit distance are in the NOrbitDistance folder.

**Run the following steps in the Windows Powershell or Linux Bash shell:**

**Step 1. Enumerate all N-Orbits in each sample and bootstrap representative vectors.**

This step generates a compiled CSV called "norbits.csv" of all N-Orbits in the dataset. It also generates a folder called "sampled_vectors" of all sampled vectors from the dataset, organized by neighborhood.

```bash
conda activate NOrbit
cd NOrbitDistance
python Step1_N-Orbit-Enumerate.py     # replace Image1 with your image name
```

**Hyperparameters**

* input_file_path: The path to your input dataset CSV file.

* intermediate_path: The path to the folder where your intermediate files will be stored. This path should be different for each dataset / analysis performed to avoid file conflicts and ensure each is performed independently.

* MODE: Set to "Neighborhood" in most cases. If neighborhood instances like with the CODEX Spleen example are desired, set to "Instance".

* im_label: The name of your image/sample label column

* x_label: The name of your x-coordinate column

* y_label: The name of your y-coordinate column

* neighborhood_label: The name of your neighborhood label column

* nucleusPenalty: Your desired nucleus change penalty *p* (recommended 1)

* radius: The maximum search radius (in microns) for cell distances to be recorded (recommended 100). Distances above this value map to a proximity score of 0.

* distMin: The minimum distance considered when measuring cell proximity (recommended 10). Everything at or below this value maps to a proximity score of 1.

* instanceRadius: Applies to "Instance" mode only. The maximum distance for two cells to be considered neighbors in radius-based spatial graph construction for instance splitting.

* minSize: The minimum number of cells in a neighborhood (within a sample) for consideration.

* sample_size: The number of bootstrapped N-Orbits used to represent each neighborhood (recommended at least 1000)

*For calculating sample-level distances, specify the neighborhood_label to the column of constant values, as mentioned in Preparing Inputs.*

**Output Files**

* norbits.csv: CSV of all N-Orbit vectors in dataset. Column labels are N-Orbit vector indices (0 through 2 x # of cell types - 1) and unit label (e.g. Image_Neighborhood).
The first half of indices correspond to the nucleus encoding (and map to cell types in alphabetical order), and the second half corresponds to the orbit encoding (also alphabetical).

* sampled_vectors: Folder of CSVs for sampled vectors of each unit. Column labels are N-Orbit vector indices.

* instance_metadata.csv (Instance Mode only): CSV with spatial coordinates and types for all cells, with an updated "unit" column for instance assignments.

* Step1_RunTime.txt: Stored runtime of Step 1.

This step takes about 20-30 minutes on the provided 100-sample synthetic dataset without parallelization.

**Hyperparameters**

* intermediate_path: The path to the folder where your intermediate files are stored, as in Step1.

**Step 2a: Compute pairwise neighborhood distances**

This step creates a "neighborhood_dists" folder, and computes a distance matrix between representative N-Orbit vectors of each neighborhood pair before calculating the minimum total distance via cost matrix optimization. This step may be parallelized in chunks by image/sample or by neighborhood, or run serially on the entire dataset **("All Mode")**. When parallelizing in neighborhood chunks **("Neighborhood Mode")**, distances are calculated from the specified neighborhood to all other neighborhoods lexicographically before it. For image/sample chunks **("Image Mode")**, the same is done for all neighborhoods in the specific image/sample. The mode is set via the UNIT_MODE parameter in the Step2a script.

For **"Neighborhood Mode"**, this step needs to be run for every neighborhood of every image in the dataset. Both the image and neighborhood labels need to be passed as a single command line argument. For example, for Neighborhood1 of Image1, the command would be as follows. 

```bash
python Step2a_Neighborhood-Distances.py Image1_Neighborhood1
```

For **"Image Mode"**, this step needs to be run for every image in the dataset. The image is passed as a command line argument as follows.

```bash
python Step2a_Neighborhood-Distances.py Image1
```

For **"All Mode"**, this step only needs to be run once, with command line argument "all".

```bash
python Step2a_Neighborhood-Distances.py all
```

**Hyperparameters**

* UNIT_MODE: The parallelization mode as just described. This can take on the values "Image", "Neighborhood", or "All".

* intermediate_path: The path to the folder where your intermediate files are stored, as earlier.

*For calculating sample-level distances, use Image Mode, substituting the constant value for the neighborhood label, e.g. Image1_0.*

This step takes a few hours on the provided 300-neighborhood synthetic dataset using bootstrap sample size 1000 and without parallelization. Runtime scales roughly cubically with N-Orbit bootstrap sample size and quadratically with neighborhood count. Parallelization is recommended for a substantially greater number of neighborhoods (1000+) and/or larger N-Orbit bootstrap sample sizes.

**Output Files**

* dists: Folder of CSVs containing column subsets of the neighborhood distance matrix for all neighborhoods specified in the command line argument, which are stored in column labels.

* RunTimeLogs: Folder of text files containing runtimes (in seconds) for each call of Step2.

**Step 2b: Compile individual distance calculation runs into a single distance matrix**

This step creates a CSV of the compiled neighborhood distance matrix from the individual CSVs in the neighborhood_dists folder. This step needs to be run regardless of UNIT_MODE used in Step2a.

```bash
python Step2b_Compile-Distance-Matrices.py
```

**Hyperparameters**

* input_file_path: The file path to the original input file, as in Step1.

* intermediate_path: The folder path where intermediates are stored, as earlier.

* im_label: The name of the image label column, as in Step1.

* neighborhood_label: The name of the neighborhood label column, as in Step1.

**Outputs**

* neighborhood_distance_matrix.csv: Final neighborhood distance matrix with columns and row labels for each unit (e.g. Image_Neigborhood).

### **N-Orbit Enrichment**

Given a neighborhood, or neighborhood cluster of interest, enriched N-Orbit may be identified using a bootstrap-permutation test. Scripts for this procedure are located in the NOrbitEnrichment folder.

**Run the following steps in the Windows Powershell or Linux Bash shell:**

**Setup**

The following should compile the Cython code and create permutation_test.c and permutation_test.html folders.

```bash
cd ../NOrbitEnrichment
python setup.py build_ext --inplace
```

**Preparing inputs**

To prepare inputs for the neighborhoods of interest, compile the corresponding CSVs in the neighborhood_vectors folder into a single CSV with an additional column annotating which image and neighborhood each vector initially came from. An example is provided in "examples/NR_cluster.csv".

**Calculating N-Orbit Enrichment**

Bootstrap-permutation tests may be run as multiple jobs to allow for parallelization. The number of tests should be kept consistent for all runs. Each run (or batch of tests) is designated through a command line argument as follows.

```bash
python n-orbit-enrichment.py 1  # example for Run 1
```

**Hyperparameters**

* vectors_path: The path to your input file storing the neighborhood vectors of interest.

* trials_path: The folder path where results will be stored for each set of trials for this neighborhood cluster.

* num_permutations: The number of bootstrap-permutation tests for each run.

* numCellTypes: The number of unique cell types in your dataset. It should be half the length of each N-Orbit vector.

* threshold: The threshold for binarization of N-Orbit vectors. Default is 0.2 which corresponds to 50 microns.

**Outputs**
A folder with the name specified in output_path with CSVs containing p and adjusted p values for each binarized vector in the set of interest.
Numerical column labels correspond to indices of the N-Orbit vector. 
Note that the "q value" column is the FDR-corrected p values for this particular run. They will be recalculated in the compile script for multiple runs.

**Compiling bootstrap permutation tests across multiple runs**

After completing all runs, the following command will compile records from trials_path into final enrichment results.

```bash
python compile.py
```

**Hyperparameters**

* input_file_path: The path to your original input files for calculating N-Orbit distance, as in Step1.

* trials_path: The folder path where results are stored for each set of trials for this neighborhood cluster, as earlier.

* cell_type_label: The name of the column for cell type labels, as in Step1.

* output_path: The file path where the final enrichment results will be stored.

**Outputs**

* n-orbit-enrichment-results.csv: A CSV of compiled enrichment results. Numerical column indices correspond to binarized N-Orbit vector indices.
The "q values" column corresponds to FDR-corrected p values calculated from the "p values" column.

## Maintainers

Barbara Xiong (barbara.xiong@pennmedicine.upenn.edu)

Yuxuan Hu (huyuxuan@xidian.edu.cn)

Kai Tan (tank1@chop.edu)

## Citation
TBD
