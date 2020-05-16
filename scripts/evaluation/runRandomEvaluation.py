import pandas as pd
import random
import sys
from rdkit.Chem import AllChem,rdReducedGraphs
import rdkit.Chem as Chem

from rdkit import DataStructs
import numpy as np
from math import sqrt
import time
from statistics import mode
import statistics
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')
import multiprocessing
import logging


def calc_ergfp( fp1, fp2 ):
    """Calulate reduced graph based similarity
    #from https://iwatobipen.wordpress.com/2016/01/16/ergfingerprint-in-rdkit/
    """
    denominator = np.sum( np.dot(fp1,fp1) ) + np.sum( np.dot(fp2,fp2) ) - np.sum( np.dot(fp1,fp2 ))
    numerator = np.sum( np.dot(fp1,fp2) )
    return numerator / denominator

def getMeasures(d):
    """clunky integration of performance measures for the validation
    @d: dictionary of true/false predictions
    return: dictionary of performance measures"""
    from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef



    tp=d["{}_ok".format(firsttype)]
    tn=d["{}_ok".format(secondtype)]
    fp=d["{}_false".format(firsttype)]
    fn=d["{}_false".format(secondtype)]

    y_true=np.hstack((np.repeat(1,tp),np.repeat(0,tn),np.repeat(0,fp),np.repeat(1,fn)))
    y_predict=np.hstack((np.repeat(1,tp),np.repeat(0,tn),np.repeat(1,fp),np.repeat(0,fn)))
    ba=balanced_accuracy_score(y_true,y_predict)
    f1=f1_score(y_true,y_predict)
    mcc=matthews_corrcoef(y_true,y_predict)
    
    return({"ba":ba,"f1":f1,"mcc":mcc,"notFound":d["notFound"],"found":d["found"]})

def getTrainTestSet(fileName,fraction,random_seed=0):
    """split training and test set
    @ filename: Filename to read smiles from
    @ fraction: fraction of dataset to build the training set, the rest will be used for testing
    @ random_seed: Optional random seed
    return: list of 2 pandas data frames, first training set, second test set
    """
    data=pd.read_csv(fileName)
    data_copy = data.copy()
    mtran=np.random.RandomState(seed=random_seed)

    train_set = data_copy.sample(frac=fraction, random_state=mtran)
    test_set = data_copy.drop(train_set.index)
    return(train_set,test_set)


def getRandomBindingMode(smiles_train):
    return(random.choice(smiles_train))
    
def evaluate(fingerprintMethod,similarityMethod,similarityThreshold=0.5,debug=0,morganFpRadius=2,evaluationMode="top"):
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
    trainDict={"type1":type1_train,"type2":type2_train,"type1_2":type1_2_train}
    testDict={"type1":type1_test,"type2":type2_test,"type1_2":type1_2_test}
    
    trainFirstTypeSmiles=trainDict["{}".format(firsttype)].smiles
    trainSecondTypeSmiles=trainDict["{}".format(secondtype)].smiles
    testFirstTypeSmiles=testDict["{}".format(firsttype)].smiles
    testSecondTypeSmiles=testDict["{}".format(secondtype)].smiles
    
    smiles_train=list(trainFirstTypeSmiles)+list(trainSecondTypeSmiles)
    smiles_test=list(testFirstTypeSmiles)+list(testSecondTypeSmiles)

    n_train_type1=len(type1_train)
    n_train_type2=len(type2_train)
    n_train_type1_2=len(type1_2_train)
    
    (fp_train,mol_train)=generateFingerprints(smiles_train,fingerprintMethod,morganFpRadius=morganFpRadius)
    (fp_test,mol_test)=generateFingerprints(smiles_test,fingerprintMethod,morganFpRadius=morganFpRadius)
    
    resList=[]
    nFound=0
    nCorrectTotal=0
    nNotFound=0
    confusionDict={"{}_ok".format(firsttype):0,"{}_false".format(firsttype):0,"{}_ok".format(secondtype):0,"{}_false".format(secondtype):0}

    #now test every molecule in the test set
    for idx,molSmile in enumerate(smiles_test):
        if fp_test[idx]!="":
            realBindingMode=getBindingModeFromSmiles(molSmile)  #get the known real binding mode of the test molecule
            randomSmiles=getRandomBindingMode(smiles_train)
            predictedBindingMode=getBindingModeFromSmiles(randomSmiles)
                    
            if predictedBindingMode==realBindingMode:
                nCorrectTotal+=1
                confusionDict[realBindingMode+"_ok"]+=1
            else :
                confusionDict[realBindingMode+"_false"]+=1

            nFound+=1
            
    confusionDict["notFound"]=nNotFound
    confusionDict["found"]=nFound
    return(confusionDict)


def generateFingerprints(smiles,fingerprintMethod,morganFpRadius=2):
    """generic fingerprint generation method that works with Morgan Fingerprints & reduced graphs
    @smiles: input smiles string
    @fingerprintMethod: pointer to fingerprint method to use to encode the molecules
    @morganFpRadius: Optional radius for morgan fingerprints
    return: list of fingerprints and corresponding rdkit molecules"""
    fp=[]
    mols=[]
    for smi in smiles:
        try:
            mol=Chem.MolFromSmiles(smi)
            if fingerprintMethod==AllChem.GetMorganFingerprint:
                fp.append(fingerprintMethod(mol,morganFpRadius))
            elif fingerprintMethod== rdReducedGraphs.GetErGFingerprint: 
                fp.append(fingerprintMethod(mol))
            mols.append(mol)
            #trainFps.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi)))
        except Exception as err:
            fp.append("")
            #print(err)
            pass
    return(fp,mols)



def getEvaluationStats(similarityMethod,similarityThreshold,n=2,morganFpRadius=2,evaluationMode="top"):
    results=[]
    for i in range(n):   
            results.append(evaluate(AllChem.GetMorganFingerprint,similarityMethod,similarityThreshold=similarityThreshold,morganFpRadius=morganFpRadius,evaluationMode=evaluationMode))

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
    print(ba,f1,mcc)
    print("ba:\t%.5f\t%.5f\tf1:\t%.5f\t%.5f\tmcc\t%.5f\t%.5f\tfound:\t%d\tnotFound\t%d"%(np.mean(ba),np.std(ba),np.mean(f1),np.std(f1),np.mean(mcc),np.std(mcc),np.mean(found),np.mean(notFound)))

if __name__ == "__main__":

    if len(sys.argv)<3:
        sys.exit("USAGE: python runRandomEvaluation.py type1 type2 [randomSeed]\nAllowed types: type1,type2,type1_2")

    firsttype: str=sys.argv[1]
    secondtype: str=sys.argv[2]
    randomSeed: int= 3333
    if(len(sys.argv)==4):
        randomSeed: int=int(sys.argv[3])
    mtran=np.random.RandomState(seed=randomSeed)

    getEvaluationStats(None,None,10,2,None)
