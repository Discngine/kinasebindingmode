import pandas as pd
import random
import sys
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit import DataStructs

def getBindingModeFromSmiles(smiles):
    if(sum(type1.smiles==smiles)): return "type1"
    if(sum(type2.smiles==smiles)): return "type2"
    if(sum(type1_2.smiles==smiles)): return "type1_2"
    

type1=pd.read_csv("prepared_data/type1.csv")
type1["type"]="type1"
type2=pd.read_csv("prepared_data/type2.csv")
type2["type"]="type2"
type1_2=pd.read_csv("prepared_data/type1_2.csv")
type1_2["type"]="type1_2"
smiles=list(type1.smiles)+list(type2.smiles)+list(type1_2.smiles)

random.shuffle(smiles)
split=int(len(smiles)/2)
train_data=smiles[:split]
trainFps=[]
for smi in train_data:
    try:
        trainFps.append(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),2))
    except Exception:
        pass

test_data=smiles[split:]
nFound=0
nNotFound=0
for idx,molSmile in enumerate(test_data):
    try:
        realBindingMode=getBindingModeFromSmiles(molSmile)
        matchingTrainIdx=[]
        matchingSimilarities=[]
        m1 = Chem.MolFromSmiles(molSmile)
        fp1 = AllChem.GetMorganFingerprint(m1,2)
        for fpidx,trainFp in enumerate(trainFps):
            similarity=DataStructs.DiceSimilarity(trainFp,fp1)
            if(similarity>0.5):
                matchingTrainIdx.append(fpidx)
                matchingSimilarities.append(similarity)
        if(len(matchingTrainIdx)):
            print("Found matching molecules for query molecule")
            print(matchingSimilarities)
            argmax=np.argmax(matchingSimilarities)
            smilesIndex=matchingTrainIdx[argmax]
            matchingSmiles=trainFps[smilesIndex]
            bindingMode=getBindingModeFromSmiles(matchingSmiles)
            print(bindingMode,realBindingMode)
            nFound+=1
        else:
            print("no match found")
            nNotFound+=1
    except Exception as err:
        pass
print(nFound)   
print(nNotFound)
    #except Exception:
    #   pass

    