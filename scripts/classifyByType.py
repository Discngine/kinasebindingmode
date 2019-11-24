import pandas as pd
import numpy as np

data=pd.read_csv("kinaseBindingModsKlifs.csv",sep="\t")

#keep only human data (for whatever non-obvious reason) - simply comment to avoid doing that.
#here I just do that to reproduce results from the Miljkovic paper
data=data.loc[(data["species"]=="Human"),:]


type1=data.loc[(data["DFG"]=="in") & data["allosteric_ligand"]=="0") & (data["aC_helix"]=="in") & pd.notna(data["smiles"]),:]
type1_2=data.loc[(data["DFG"]=="in") & data["allosteric_ligand"]=="0") & (data["aC_helix"]=="out") & pd.notna(data["smiles"]),:]
type2=data.loc[(data["DFG"]=="out") & data["allosteric_ligand"]=="0") & (data["aC_helix"]=="out") & (data["back"]) & pd.notna(data["smiles"]),:]
allosteric=data.loc[(data["allosteric_ligand"]!="0"),:] #issue with 4oli -> revalidate / classify binding sites here before using this

#number of structures with NON allosteric inhibitors
print(len(data.loc[(data["allosteric_ligand"]=="0")  & pd.notna(data["smiles"]),"pdb"].unique()))


#getting unique ligands for different types of binders
uniqueType1Smiles=type1.smiles.unique() #1941 unique smiles 30th October 2019
print(len(uniqueType1Smiles),"unique type 1 molecules")
uniqueType1_2Smiles=type1_2.smiles.unique() #592 unique smiles 30th October 2019
print(len(uniqueType1_2Smiles), "unique type 1 1/2 molecules")
uniqueType2Smiles=type2.smiles.unique() #270 unique smiles 30th October 2019
print(len(uniqueType2Smiles),"unique type 2 molecules")

tmp=[ligand for ligand in type1.ligand if sum( type2.ligand==ligand) or sum(type1_2.ligand==ligand) ]
tmp+=[ligand for ligand in type2.ligand if sum( type1.ligand==ligand) or sum(type1_2.ligand==ligand) ]
tmp+=[ligand for ligand in type1_2.ligand if sum( type1.ligand==ligand) or sum(type2.ligand==ligand) ]
ligandsToDiscard=np.unique(tmp)

#let's get ligands where the binding mode is unclear (common compounds between lists)

for ligand in ligandsToDiscard:
    type1 = type1[type1.ligand != ligand]
    type1_2 = type1_2[type1_2.ligand != ligand]
    type2 = type2[type2.ligand != ligand]


#getting unique ligands for different types of binders
uniqueType1Smiles=type1.smiles.unique() #1941 unique smiles 30th October 2019
print(len(uniqueType1Smiles),"unique type 1 molecules")
uniqueType1_2Smiles=type1_2.smiles.unique() #592 unique smiles 30th October 2019
print(len(uniqueType1_2Smiles), "unique type 1 1/2 molecules")
uniqueType2Smiles=type2.smiles.unique() #270 unique smiles 30th October 2019
print(len(uniqueType2Smiles),"unique type 2 molecules")

unique_type1=type1.drop_duplicates(subset="smiles")
unique_type2=type2.drop_duplicates(subset="smiles")
unique_type1_2=type1_2.drop_duplicates(subset="smiles")


unique_type1.to_csv("prepared_data/type1.csv")
unique_type1_2.to_csv("prepared_data/type1_2.csv")
unique_type2.to_csv("prepared_data/type2.csv")
#try to set up global model dataset approach according to Miljkovic et al:



#--> I cannot reproduce numbers from Miljkovic et al here. Don't understand how they 
