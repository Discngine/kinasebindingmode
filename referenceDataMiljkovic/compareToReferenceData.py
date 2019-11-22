import pandas as pd 
import numpy as np
refdata=pd.read_table("data.dat")
ourdata=pd.read_table("../kinaseBindingModesKlifs.csv")


#first check if all of the reference molecules are found in our data
for idx,rescode in enumerate(ourdata["ligand"]):
    if rescode==0:
        rescode=ourdata.iloc[idx]["allosteric_ligand"]
    if rescode!=0 and rescode !="0":
        match=refdata["cid"].str.find(rescode.upper())
        if len(match.loc[match>-1])>0:
            print(rescode+' found')
        else:
            print(rescode+" not found")
