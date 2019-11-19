# README

## How to install the development environment
Here I used conda (miniconda, anaconda, whatever you prefer
Simply install the conda environment on your machine using.
`conda env create -f kinaseBindingMode.yml`

## Activate the environment
`conda activate kinaseBindingMode

## Data
- kinaseBindingModesKlifs.csv contains a raw extraction from klifs (https://klifs.vu-compmedchem.nl/) and the RCSB PDB (https://www.rcsb.org/). It's a plain extraction from KLIFS without any modification.
- establishReferenceSet.py is the python script to be used to build the aforementioned csv file.
