import glob, csv, sys, time
import pandas as pd
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import urllib3
import urllib.parse
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import numpy as np
from rdkit import Chem




def getLigandStructure(structure_id,session):
    """Retrieve mol2 structure of the ligand for one particular structure id and 
    transform this to a smiles string"""
    query="https://klifs.vu-compmedchem.nl/api/structure_get_ligand?structure_ID="+str(structure_id)
    response=session.get(query)
    if response.status_code==200:
        try:
            mol=Chem.rdmolfiles.MolFromMol2Block(response.content)
            return(Chem.MolToSmiles(mol))
        except Exception as err:
            print(err)
            return("NA")



#let's get all structures available in KLIFS currently

sklifs = Session()
klifsHost="http://klifs.vu-compmedchem.nl/api/"

#get all kinase names :
endpoint=klifsHost+"/kinase_names"
response=sklifs.get(endpoint)

#prepare outputfile handle so we can reuse what we extract here later for training etc
of=open("kinaseBindingModsKlifs.csv","w")
of.write("Kinase ID\tKinase Name\tsmiles\tstructure_ID\tpdb\talt\tchain\tmissing_residues\tligand\tallosteric_ligand\tDFG\taC_helix\tback\n")
if response.status_code==200:
    for kinase in response.json():
        kinaseID=kinase["kinase_ID"]
        species=kinase["species"]
        print("Analyzing Kinase "+kinase["name"])
        #if kinaseID==441:
        endpoint=klifsHost+"/structures_list?kinase_ID="+str(kinaseID)
        structureResponse=sklifs.get(endpoint)

        descriptors=["structure_ID","pdb","alt","chain","missing_residues","ligand","allosteric_ligand","DFG","aC_helix","back"]
        if structureResponse.status_code==200:
            
            structureList=structureResponse.json()
            for structure in structureList:
                structureID=structure["structure_ID"]
                smiles=getLigandStructure(structureID,sklifs)

                of.write(str(kinaseID)+"\t"+kinase["name"]+"\t"+smiles)
                for descriptor in descriptors:
                    of.write("\t"+str(structure[descriptor]))
                of.write("\n")
of.close()

