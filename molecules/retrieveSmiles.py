import pandas as pd 
import numpy as np
data=pd.read_csv("../prepared_data/kinaseBindingModesKlifs_andKeyResidues.csv",sep=";")

het=pd.read_csv("ambiguous.txt")

for code in het["residueCode"]:
    smi=list((data.loc[data['ligand'] == code,"smiles"]))
    print(smi[0])
    het.loc[het["residueCode"]==code,"smiles"]=smi[0]

het.to_csv("ambiguous.csv",index=False)