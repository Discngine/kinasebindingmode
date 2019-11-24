# README

## How to install the development environment
Here I used conda (miniconda, anaconda, whatever you prefer).
Simply install the conda environment on your machine using.
`conda env create -f kinaseBindingMode.yml`

## Activate the environment
`conda activate kinaseBindingMode`

## Scripts:
`evaluateGlobalSets.py:` contains the 50/50 evaluation protocol

## Folder structure
### referenceDataMiljkovic
Contains the reference data from Miljkovic et al 2019 in data.dat. It also contains a python script comparing the data used here vs theirs. 
### prepared_data
Contains prepared smiles strings separated by binding Mode (type1, type2, type1_2). 
kinaseBindingModesKlifs_andKeyResidues.csv contains a raw extraction from KLIFS, the PDB and 3decision (for key residue identification)


## Data
- kinaseBindingModesKlifs.csv contains a raw extraction from klifs (https://klifs.vu-compmedchem.nl/) and the RCSB PDB (https://www.rcsb.org/). It's a plain extraction from KLIFS without any modification.
- establishReferenceSet.py is the python script to be used to build the aforementioned csv file.
