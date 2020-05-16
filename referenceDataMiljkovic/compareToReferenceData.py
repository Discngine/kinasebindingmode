import pandas as pd 
import numpy as np
refdata=pd.read_table("data.dat")
ourdata=pd.read_table("../scripts/kinaseBindingModesKlifs.csv")

if 0 : 

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


#check if in the reference data some stuctures are in there and not in the Miljkovic data? 
found=0
notfound=0
for idx,cid in enumerate(refdata["cid"]):
    if(len(cid.split(";"))>1):
        cids=cid.split(";")
        for cidsplit in cids:
            cidsplit=cidsplit.strip()
            match=ourdata["ligand"].str.find(cidsplit.upper())
            match2=ourdata["allosteric_ligand"].str.find(cidsplit.upper())
            if len(match.loc[match>-1]) or len(match.loc[match2>-1]):
                found+=1
            else :
                notfound+=1
                print(cid)
    else:
        match=ourdata["ligand"].str.find(cid.upper())
        match2=ourdata["allosteric_ligand"].str.find(cid.upper())
        if len(match.loc[match>-1]) or len(match.loc[match2>-1]):
            found+=1
        else :
            notfound+=1
            print(cid)
print(found)
print(notfound)