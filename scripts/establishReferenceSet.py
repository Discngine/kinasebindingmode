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
import xmltodict
import time


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



def getReleaseDateFromPDB(pdbCode,session):
    """Retrieve mol2 structure of the ligand for one particular structure id and 
    transform this to a smiles string"""
    query="https://www.rcsb.org/pdb/rest/getEntityInfo?structureId="+str(pdbCode)
    response=session.get(query)
    if response.status_code==200:
        try:
            d=xmltodict.parse(response.content)
            return(d["entityInfo"]["PDB"]["@release_date"])
        except Exception as err:
            print("Error fetching release date from PDB for pdb code: "+str(pdbCode))
            print(err)
            return("NA")

def getLigandStructureFromPDB(het,session):
    """Retrieve mol2 structure of the ligand for one particular structure id and 
    transform this to a smiles string"""
    query="https://www.rcsb.org/pdb/rest/describeHet?chemicalID="+str(het)
    response=session.get(query)
    if response.status_code==200:
        try:
            d=xmltodict.parse(response.content)
            return(d["describeHet"]["ligandInfo"]["ligand"]["smiles"])
        except Exception as err:
            print("Error fetching smiles from PDB for het code: "+str(het))
            print(err)
            return("NA")


#let's get all structures available in KLIFS currently

sklifs = Session()
spdb = Session()
klifsHost="http://klifs.vu-compmedchem.nl/api/"

#get all kinase names :
endpoint=klifsHost+"/kinase_names"
response=sklifs.get(endpoint)

#just not to call the PDB unnecessarily
ligandDict={}


#prepare outputfile handle so we can reuse what we extract here later for training etc
of=open("kinaseBindingModesKlifs.csv","w")
of.write("Kinase ID\tKinase Name\tsmiles\tstructure_ID\tpdb\talt\tchain\tmissing_residues\tligand\tallosteric_ligand\tDFG\taC_helix\tback\tspecies\treleaseDate\n")
if response.status_code==200:
    for kinase in response.json():
        kinaseID=kinase["kinase_ID"]
        species=kinase["species"]
        print("Analyzing Kinase "+kinase["name"])
        #if kinaseID==441:
        endpoint=klifsHost+"/structures_list?kinase_ID="+str(kinaseID)
        try:
            structureResponse=sklifs.get(endpoint)
        except Exception as err:
            print(err)
            time.sleep(10)
            structureResponse=sklifs.get(endpoint)
            pass

        descriptors=["structure_ID","pdb","alt","chain","missing_residues","ligand","allosteric_ligand","DFG","aC_helix","back"]
        if structureResponse.status_code==200:
            
            structureList=structureResponse.json()
            for structure in structureList:
                structureID=structure["structure_ID"]
                #smiles=getLigandStructure(structureID,sklifs)
                if len(str(structure["ligand"]))==3 and structure["ligand"]: 
                    if structure["ligand"] not in ligandDict:
                        smiles=getLigandStructureFromPDB(structure["ligand"],spdb)
                        ligandDict[structure["ligand"]]=smiles
                    else :
                        smiles=ligandDict[structure["ligand"]]
                else :
                    smiles="NA"
                if smiles is None:
                    smiles="NA"
                releaseDate=getReleaseDateFromPDB(structure["pdb"],spdb)
                if releaseDate is None:
                    releaseDate="NA"
                print(kinase["name"])
                print(smiles)
                of.write(str(kinaseID)+"\t"+kinase["name"]+"\t"+smiles)
                for descriptor in descriptors:
                    of.write("\t"+str(structure[descriptor]))
                of.write("\t"+species+"\t"+releaseDate+"\n")


of.close()

