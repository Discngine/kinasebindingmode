import pandas as pd
import random
import sys
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit import DataStructs
import numpy as np
from rdkit.Chem import MACCSkeys

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
smiles=list(type1.smiles)+list(type2.smiles)#+list(type1_2.smiles)

random.shuffle(smiles)
split=int(len(smiles)/1.7)
train_data=smiles[:split]
trainFps=[]
for smi in train_data:
    try:
        trainFps.append(AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),3))
        #trainFps.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi)))
    except Exception:
        trainFps.append("")
        pass

test_data=smiles[split:]
nFound=0
nCorrectTotal=0
nNotFound=0
confusionDict={"type1_tp":0,"type1_fn":0,"type2_tp":0,"type2_fn":0}
for idx,molSmile in enumerate(test_data):
    try:
        realBindingMode=getBindingModeFromSmiles(molSmile)
        
        matchingTrainIdx=[]
        matchingSimilarities=[]
        m1 = Chem.MolFromSmiles(molSmile)
        fp1 = AllChem.GetMorganFingerprint(m1,3)
        #fp1 = MACCSkeys.GenMACCSKeys(m1)
        for fpidx,trainFp in enumerate(trainFps):
            if(trainFp!=""):
                #similarity=DataStructs.FingerprintSimilarity(trainFp,fp1)
                similarity=DataStructs.DiceSimilarity(trainFp,fp1)
                if(similarity>0.5):
                    matchingTrainIdx.append(fpidx)
                    matchingSimilarities.append(similarity)
                    #if realBindingMode=="type2":
                    #    print(molSmile,train_data[fpidx])
        if(len(matchingTrainIdx)):
            nCorrect=0
            for matchingTrainIndex in matchingTrainIdx:
            
                matchingSmiles=train_data[matchingTrainIndex]
                bindingMode=getBindingModeFromSmiles(matchingSmiles)
                #print(bindingMode,realBindingMode)
                if bindingMode==realBindingMode:
                    nCorrect+=1
                    
            #print(nCorrect/float(len(matchingTrainIdx)))
            if(nCorrect/float((len(matchingTrainIdx)))>0.5):
                nCorrectTotal+=1
                confusionDict[realBindingMode+"_tp"]+=1
            else:
                confusionDict[realBindingMode+"_fn"]+=1
                    


            nFound+=1
        else:
            
            nNotFound+=1
    except Exception as err:
        pass
print("Found")
print(nFound)
print("correct") 
print(nCorrectTotal)
print("not found")
print(nNotFound)
print(confusionDict)
    #except Exception:
    #   pass

    