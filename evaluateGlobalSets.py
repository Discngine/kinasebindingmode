import pandas as pd
import random
import sys
from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from rdkit.Chem.Fraggle.FraggleSim import GetFraggleSimilarity
from rdkit import DataStructs
import numpy as np
from rdkit.Chem import MACCSkeys
from math import sqrt
import time
from statistics import mode
import statistics
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

mtran=np.random.RandomState(seed=3333)

def getMeasures(d):
    tp=d["type1_ok"]
    tn=d["type2_ok"]
    fp=d["type1_false"]
    fn=d["type2_false"]
    if((d["type1_ok"]+d["type1_false"]))>0.0:
        ppv=d["type1_ok"]/(d["type1_ok"]+d["type1_false"])
    else:
        ppv=0.0
    if (d["type1_ok"]+d["type2_false"])>0.0:
        tpr=d["type1_ok"]/(d["type1_ok"]+d["type2_false"])
    else:
        tpr=0.0
    if (d["type2_ok"]+d["type1_false"])>0.0:
        tnr=d["type2_ok"]/(d["type2_ok"]+d["type1_false"])
    else:
        tnr=0.0
    ba=(tpr+tnr)/2.0
    if (ppv+tpr)>0.0:
        f1=(2.0*ppv*tpr)/(ppv+tpr)
    else :
        f1=0.0
    if((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)>0.0):
        mcc=(tp*tn-fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    else :
        mcc=0.0
    return({"ba":ba,"f1":f1,"mcc":mcc,"notFound":d["notFound"],"found":d["found"]})

def getTrainTestSet(fileName,fraction,random_seed=0):
    data=pd.read_csv(fileName)
    data_copy = data.copy()
    train_set = data_copy.sample(frac=fraction, random_state=mtran)
    test_set = data_copy.drop(train_set.index)
    return(train_set,test_set)

"""TODO: introduce pointers to functions for fingerprints"""
def generateFingerprints(smiles,fingerprintMethod):
    fp=[]
    for smi in smiles:
        try:
            fp.append(fingerprintMethod(Chem.MolFromSmiles(smi),2))
            #trainFps.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi)))
        except Exception as err:
            fp.append("")
            pass
    return(fp)


def evaluate(fingerprintMethod,similarityMethod,similarityThreshold=0.5,acceptanceRatio=0.5,debug=0):
    def getBindingModeFromSmiles(smiles):
        if(sum(type1_train.smiles==smiles) or sum(type1_test.smiles==smiles)): return "type1"
        if(sum(type2_train.smiles==smiles) or sum(type2_test.smiles==smiles)): return "type2"
        if(sum(type1_2_train.smiles==smiles) or sum(type1_2_test.smiles==smiles)): return "type1_2"



    (type1_train,type1_test)=getTrainTestSet("prepared_data/type1.csv",0.5)
    type1_train["type"]="type1"
    type1_test["type"]="type1"
    if debug:    
        print("TYPE 1:")
        print("      - %d molecules in training set."%(len(type1_train)))
        print("      - %d molecules in test set."%(len(type1_test)))
        
    
    (type2_train,type2_test)=getTrainTestSet("prepared_data/type2.csv",0.5)
    type2_train["type"]="type2"
    type2_test["type"]="type2"

    if debug:  
        print("TYPE 2:")
        print("      - %d molecules in training set."%(len(type2_train)))
        print("      - %d molecules in test set."%(len(type2_test)))
        
    
    (type1_2_train,type1_2_test)=getTrainTestSet("prepared_data/type1_2.csv",0.5)
    type1_2_train["type"]="type1_2"
    type1_2_test["type"]="type1_2"
    
    if debug:  
        print("TYPE 1 1/2:")
        print("      - %d molecules in training set."%(len(type1_2_train)))
        print("      - %d molecules in test set."%(len(type1_2_test)))
        
    
    smiles_train=list(type1_train.smiles)+list(type2_train.smiles)
    smiles_test=list(type1_test.smiles)+list(type2_test.smiles)

    n_train_type1=len(type1_train)
    n_train_type2=len(type2_train)
    n_train_type1_2=len(type1_2_train)
    
    fp_train=generateFingerprints(smiles_train,fingerprintMethod)
    fp_test=generateFingerprints(smiles_test,fingerprintMethod)
    
    resList=[]
    nFound=0
    nCorrectTotal=0
    nNotFound=0
    confusionDict={"type1_ok":0,"type1_false":0,"type2_ok":0,"type2_false":0}

    #now test every molecule in the test set
    for idx,molSmile in enumerate(smiles_test):
        #only if the fingerprint could be calculated
        if fp_test[idx]!="":
            realBindingMode=getBindingModeFromSmiles(molSmile)  #get the known real binding mode of the test molecule
            
            matchingTrainIdx=[] #list to contain all matching molecules with a similarity above a given threshold
            matchingSimilarities=[] #list to get actual similarities
            fp1=fp_test[idx]
            
            #search if there's a similar molecule in the training set or not. If yes, check its binding mode
            for fpidx,trainFp in enumerate(fp_train):
                if(trainFp!=""):
                    #similarity=DataStructs.FingerprintSimilarity(trainFp,fp1)
                    
                    if similarityMethod.__name__=="GetFraggleSimilarity":
                        try:
                            similarity,match=similarityMethod(Chem.MolFromSmiles(molSmile),Chem.MolFromSmiles(smiles_train[fpidx]))
                        except Exception:
                            print("Failed for smiles : "+smiles_train[fpidx]+" and "+molSmile)
                            similarity=0.0
                            pass

                    else:
                        similarity=similarityMethod(trainFp,fp1)
                    
                    if(similarity>=similarityThreshold):
                        matchingTrainIdx.append(fpidx)
                        matchingSimilarities.append(similarity)

            #if at least some matches were found
            if(len(matchingTrainIdx)):
                nCorrect=0
                #count how many times we observe the same binding mode as the known real binding mode among the matching molecules
                matchingBindingModes=[]
                for matchingTrainIndex in matchingTrainIdx:
                
                    matchingSmiles=smiles_train[matchingTrainIndex]
                    bindingMode=getBindingModeFromSmiles(matchingSmiles)
                    matchingBindingModes.append(bindingMode)
                    #print(bindingMode,realBindingMode)
                    if bindingMode==realBindingMode:
                        nCorrect+=1
                        
                
                #if more than x of these matching molecules are correct we consider that the similarity search result is "predictive"
                #this is very biased ... in a proper prediction setting, here we need to define what that actually means. 
                #Explanation: if acceptanceRatio >0.6, this means that in accepted similar molecules, more than 60% are in the correct class.
                try:
                    a = np.array(matchingBindingModes)
                    unique, counts = np.unique(a, return_counts=True)
                    bindingModeData = pd.DataFrame({'mode':unique, 'counts':counts})
                    bindingModeData.loc[(bindingModeData["mode"]=="type1"),"counts"]/=n_train_type1
                    bindingModeData.loc[(bindingModeData["mode"]=="type2"),"counts"]/=n_train_type2
                    bindingModeData.loc[(bindingModeData["mode"]=="type1_2"),"counts"]/=n_train_type1_2

                    # for bidx,bindingMode in enumerat(unique):
                    #     if bindingMode=="type1":
                            
                    #     elif bindingMode=="type2":

                    #     elif bindingMode=="type1_2":
                    
                    #matchingBindingMode=mode(matchingBindingModes)
                    predictedBindingMode=bindingModeData.loc[bindingModeData['counts'].idxmax(),"mode"]
                    if predictedBindingMode==realBindingMode:
                        nCorrectTotal+=1
                        confusionDict[realBindingMode+"_ok"]+=1
                    else :
                        confusionDict[realBindingMode+"_false"]+=1
                
                except statistics.StatisticsError:
                    confusionDict[realBindingMode+"_false"]+=1
                
                
                #if(nCorrect/float((len(matchingTrainIdx)))>=acceptanceRatio):
                #    nCorrectTotal+=1
                #    confusionDict[realBindingMode+"_ok"]+=1
                #else:
                #    confusionDict[realBindingMode+"_false"]+=1

                nFound+=1
            else:
                #if no similar ligands were found, then keep track of that
                nNotFound+=1
    confusionDict["notFound"]=nNotFound
    confusionDict["found"]=nFound
    return(confusionDict)


if __name__ == "__main__":
    similarityMethods=[GetFraggleSimilarity,DataStructs.TanimotoSimilarity,DataStructs.DiceSimilarity]
    #similarityMethods=[DataStructs.TanimotoSimilarity,DataStructs.FingerprintSimilarity,DataStructs.DiceSimilarity,DataStructs.AsymmetricSimilarity,DataStructs.BraunBlanquetSimilarity,DataStructs.CosineSimilarity,DataStructs.KulczynskiSimilarity,DataStructs.McConnaugheySimilarity, DataStructs.RogotGoldbergSimilarity, DataStructs.RusselSimilarity,DataStructs.SokalSimilarity,DataStructs.TverskySimilarity]
    #acceptanceRatios=np.arange(0.1,1.1,0.1)
    acceptanceRatios=[0.1]
    similarityThresholds=np.arange(0.9,0.1,-0.05)
    for similarityMethod in similarityMethods:
        for acceptanceRatio in acceptanceRatios:
            for similarityThreshold in similarityThresholds:
                    results=[]
                    
                #try:
                    for i in range(2):    
                        results.append(evaluate(AllChem.GetMorganFingerprint,similarityMethod,acceptanceRatio=acceptanceRatio,similarityThreshold=similarityThreshold))

                    ba=[]
                    f1=[]
                    mcc=[]
                    notFound=[]
                    found=[]

                    for result in results:
                        m=getMeasures(result)
                        ba.append(m["ba"])
                        f1.append(m["f1"])
                        mcc.append(m["mcc"])
                        notFound.append(m["notFound"])
                        found.append(m["found"])
                        
                    print("%s : similarityThreshold: %.2f  acceptanceRatio: %.2f   ba: %.5f+/-%.5f   f1: %.5f+/-%.5f   mcc: %.5f+/-%.5f   found: %d+/-%.2f    notFound %d+/-%.2f"%(similarityMethod.__name__,similarityThreshold,acceptanceRatio,np.mean(ba),np.std(ba),np.mean(f1),np.std(f1),np.mean(mcc),np.std(mcc),np.mean(found),np.std(found),np.mean(notFound),np.std(notFound)))
                #except Exception as err:
                #    print(similarityMethod.__name__+" no results")
                #    print(err)
                #    pass
