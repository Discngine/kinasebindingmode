import pandas as pd

data=pd.read_csv("kinaseBindingModsKlifs.csv",sep="\t")

type1=data.loc[(data["DFG"]=="in") & (data["aC_helix"]=="in") & pd.notna(data["smiles"]),:]
type1_2=data.loc[(data["DFG"]=="in") & (data["aC_helix"]=="out") & pd.notna(data["smiles"]),:]
type2=data.loc[(data["DFG"]=="out") & (data["aC_helix"]=="out") & (data["back"]) & pd.notna(data["smiles"]),:]
allosteric=data.loc[(data["allosteric_ligand"]!="0"),:] #issue with 4oli -> revalidate / classify binding sites here before using this


#getting unique ligands for different types of binders
uniqueType1Smiles=type1.smiles.unique()
uniqueType1_2Smiles=type1_2.smiles.unique()
uniqueType2Smiles=type2.smiles.unique()
