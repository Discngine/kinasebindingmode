# README

## How to install the development environment
Here I used conda (miniconda, anaconda, whatever you prefer).
Simply install the conda environment on your machine using.
`conda env create -f kinaseBindingMode.yml`

## Activate the environment
`conda activate kinaseBindingMode`

## Scripts:

### I: Prepare initial public data-set
Download all annotations from KLIFS and gather ligand structures from the RCSB PDB and writes them to scripts/kinaseBindingModesKlifs.csv. This is raw annotation taken as input for other scripts. You should be able to run this script at any later time to get an updated data-set.

`cd scripts; python establishReferenceSet.py`

Next rules to split all inhibitors into different inhibitor type classes are available are encoded in the scripts/classifyByType.py

`cd scripts; python classifyByType.py`

This classfication step does a bit of cleaning un duplicate checking and puts output into prepared_data/typeX.csv files. We followed here the rules used by Miljkovic et et al.
At the time of writing the paper the following statistics were obtained: 
- 1775 unique type 1 molecules
- 495 unique type 1 1/2 molecules
- 238 unique type 2 molecules

### II: Comparing input data to Miljkovich data-set
In the folder referenceDataMiljkovic there's a simple script to compare our dataset to the one from Miljkovic. At the time of the paper this scripts outputs two molecules not corresponding between their data-set and ours: 
- residue code: 75F
- residue codes: X4Z;7KC

A total of 2075 molecules is contained in their data-set vs 2508 in ours. 


### III: Running 50/50 global evaluation
Run the global evaluation usnig the 50/50 strategy. You can run this using the runEvaluation python script such as:

`python scripts/evaluation/runEvaluation.py typeX typeY [ncpus] [randomSeed]`

ncpus and randomseed are optional arguments. 
- typeX: either **type1**, **type2** or **type1_2**
- typeY: either **type1**, **type2** or **type1_2** (not the same as typeX)

The script will run evaluations using different fingerprint radii, similarity metrics etc to distinguish typeX from typeY inhibitors and print this to the stdout. Make sure to redirect the results to a file if you want to persist them.
Example to get the performance of a type1 vs type2 classifier: 

`python scripts/evaluation/runEvaluation.py type1 type2 12 3333`

### IV: Get Random Evaluation
In order to compare previous obtained predictions to random, we've provided a script that runs a random selection for all possible classifiers: 

`python scripts/evaluation/runRandomEvaluation.py type1 type2` 

will calculate the balanced accuracy, F1 and Mathews correlation coefficient for a type1 type2 random classifier with the actual data (biased). 



## Folder structure
### referenceDataMiljkovic
Contains the reference data from Miljkovic et al 2019 in data.dat. It also contains a python script comparing the data used here vs theirs. 
### prepared_data
Contains prepared smiles strings separated by binding Mode (type1, type2, type1_2). 
kinaseBindingModesKlifs_andKeyResidues.csv contains a raw extraction from KLIFS, the PDB and 3decision (for key residue identification)


## Data
- kinaseBindingModesKlifs.csv contains a raw extraction from klifs (https://klifs.vu-compmedchem.nl/) and the RCSB PDB (https://www.rcsb.org/). It's a plain extraction from KLIFS without any modification.
- establishReferenceSet.py is the python script to be used to build the aforementioned csv file.
