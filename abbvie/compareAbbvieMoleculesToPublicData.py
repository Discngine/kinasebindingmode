import pandas as pd
import numpy as np 
from rdkit.Chem import AllChem
from rdkit import Chem



def getInchiFromSmiles(input):
    inchi=[]
    for smiles in input:
        try:
            mol=AllChem.MolFromSmiles(smiles)
            
            Chem.SanitizeMol(mol)
            inchikey=AllChem.MolToInchi(mol)
            inchi.append(inchikey)
            #print(inchikey)
        except Exception :
            print("FAILED: ", smiles)
            pass
    return(inchi)


data=pd.read_csv("/Users/peter/Desktop/AbbvieKinConformations.csv",sep=";")
abv_smiles=np.unique(data.loc[~data["smiles"].isnull(),"smiles"])
abv_inchi=getInchiFromSmiles(abv_smiles)


type1=pd.read_csv("../prepared_data/type1.csv")
type1_smiles=np.unique(type1.loc[~type1["smiles"].isnull(),"smiles"])
pdb_inchi=getInchiFromSmiles(type1_smiles)

type2=pd.read_csv("../prepared_data/type2.csv")
type2_smiles=np.unique(type2.loc[~type2["smiles"].isnull(),"smiles"])
pdb_inchi+=getInchiFromSmiles(type2_smiles)

type1_2=pd.read_csv("../prepared_data/type1_2.csv")
type1_2_smiles=np.unique(type1_2.loc[~type1_2["smiles"].isnull(),"smiles"])
pdb_inchi+=getInchiFromSmiles(type1_2_smiles)

nCommon=0
for abv in abv_inchi:
    if abv in pdb_inchi:
        nCommon+=1
print(len(abv_inchi))
print(nCommon)
print(float(nCommon)/len(abv_inchi))