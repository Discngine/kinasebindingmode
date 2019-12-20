getKinaseStrAnnotations
=======================


This is a simple script that read a kinase pdb structure, and based on specific key residues provided as option, return structural annotations. 

Requirements
------------
the script run within python3 environment. Required libraries are specified in Requirement.txt

Installation
------------
with conda :
``conda install --yes --file requirements.txt``
or if you are not using conda :
``python3 -m pip install -r requirements.txt``

How to execute the script
-------------------------

to get help:
``python getKinaseProperties.py  -h``

ex of command line:
``python getKinaseProperties.py  1:B:194 1:B:195 1:B:90 1:B:74 ../../prepared_data/structures/5l4q.pdb``

use awk to apply the script on the 60 first dataset:
``awk -F '\t' ' NR > 0  { if(NR==1) print "#",$1,$5,$6,$7,$11,$12,"predicted_DFG_and_HaC"; if(NR>1){  "python ../scripts/getKinaseStrAnnotations/getKinaseStrAnnotations.py   "$21" "$22 | getline toto; print NR, $1, $5, $6, $7, $11, $12, toto } } ' kinaseBindingModesKlifs_andKeyResidues.csv > kinaseBindingModesKlifs_andKeyResidues.out``
