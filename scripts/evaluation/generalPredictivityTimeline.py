import pandas as pd
import numpy as np 
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('max_colwidth', -1)
debug=0
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')    
import time
import sklearn
import random
import sys
from rdkit import DataStructs
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


simMetric="tanimoto"

if(len(sys.argv)<3):
    sys.exit("USAGE: python generalPredictivity.py similarityMetric outputFilePrefix\nAllowed metrics: tanimoto,dice,tversky")

simMetric:str=sys.argv[1]
outputPrefix:str=sys.argv[2]


def getMaxSimilarity(test_fp,train_fp,train_outcome):
    response=[]
    similarity=[]
    for ix,fp in enumerate(test_fp):
        maxSim=0.0
        curRefIx=-1
        if(fp!=''):
            for ref_ix,ref_fp in enumerate(train_fp):
                if(ref_fp!=''):
                    if simMetric=="tanimoto":
                        sim=DataStructs.TanimotoSimilarity(fp,ref_fp)
                    if simMetric=="tversky":
                        sim=DataStructs.TverskySimilarity(fp,ref_fp,0.9,0.1)
                    if simMetric=="dice":
                        sim=DataStructs.DiceSimilarity(fp,ref_fp)
                    if(sim>maxSim):
                        curRefIx=ref_ix
                        maxSim=sim
            predictedBindingMode=train_outcome[curRefIx]
            similarity.append(maxSim)
            response.append(predictedBindingMode)
        else: 
            response.append(0)
            similarity.append(0)
    return (response,similarity)


def getRandomPrediction(test_fp,train_fp,train_outcome):
    response=[]
    similarity=[]
    for ix,fp in enumerate(test_fp):
        maxSim=0.0
        curRefIx=-1
        if(fp!=''):
            ref_fp=''
            while(ref_fp=='') :
                ref_ix=random.choice(range(len(train_fp)))
                ref_fp=train_fp[ref_ix]
            if simMetric=="tanimoto":
                sim=DataStructs.TanimotoSimilarity(fp,ref_fp)
            if simMetric=="tversky":
                sim=DataStructs.TverskySimilarity(fp,ref_fp,0.9,0.1)
            if simMetric=="dice":
                sim=DataStructs.DiceSimilarity(fp,ref_fp)
    
            predictedBindingMode=train_outcome[ref_ix]
            similarity.append(sim)
            response.append(predictedBindingMode)
        else: 
            response.append(0)
            similarity.append(0)
    return (response,similarity)



def getRandomSet(input_fp,choices):
    response=[]
    for ix,abv_fp in enumerate(input_fp):
        response.append(random.choice(choices))
    return response

def getBindingModeFromSmiles(smiles):
        if(sum(type1_train.smiles==smiles) or sum(type1_test.smiles==smiles)): return "type1"
        if(sum(type2_train.smiles==smiles) or sum(type2_test.smiles==smiles)): return "type2"
        if(sum(type1_2_train.smiles==smiles) or sum(type1_2_test.smiles==smiles)): return "type1_2"


def getTrainTestSet(fileName,currentDate,random_seed=None):
    """split training and test set
    @ filename: Filename to read smiles from
    @ fraction: fraction of dataset to build the training set, the rest will be used for testing
    @ random_seed: Optional random seed
    return: list of 2 pandas data frames, first training set, second test set
    """
    data=pd.read_csv(fileName)
    data["releaseDate"]=pd.to_datetime(data["releaseDate"],infer_datetime_format=True)
    train_set=data[data["releaseDate"]<currentDate]
    test_set=data[data["releaseDate"]>=currentDate]
    
    return(train_set,test_set)


def generateFingerprints(smiles,fingerprintMethod,morganFpRadius=2):
    """generic fingerprint generation method that works with Morgan Fingerprints & reduced graphs
    @smiles: input smiles string
    @fingerprintMethod: pointer to fingerprint method to use to encode the molecules
    @morganFpRadius: Optional radius for morgan fingerprints
    return: list of fingerprints and corresponding rdkit molecules"""
    fp=[]
    mols=[]
    sm=[]
    for smi in smiles:
        try:
            mol=Chem.MolFromSmiles(smi)
            if fingerprintMethod==AllChem.GetMorganFingerprint:
                fp.append(fingerprintMethod(mol,morganFpRadius))
            elif fingerprintMethod== rdReducedGraphs.GetErGFingerprint: 
                fp.append(fingerprintMethod(mol))
            mols.append(mol)
            sm.append(smi)
            #trainFps.append(MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi)))
        except Exception as err:
            fp.append("")
            mols.append("")
            sm.append(smi)
            #print(err)
            pass
    return(fp,mols,sm)


for radius in range(1,4):
    from datetime import date
    from dateutil.relativedelta import relativedelta
    start_date = date(1999, 1, 1)
    end_date = date(2019, 12, 31)
    delta = relativedelta(months=+6)

    o=open("{}_radius_{}_{}.csv".format(outputPrefix,radius,simMetric),"w")
    o.write("Cycle\tCorrect\tPrediction\tActualClass\tSimilarity\tsmiles\trandomCorrect\trandomSimilarity\tdateThreshold\n")


    while start_date <= end_date:
        print (start_date.strftime("%Y-%m-%d"))
        start_date += delta

        for cycle in range(10):

            (type1_train,type1_test)=getTrainTestSet("prepared_data/type1.csv",start_date)

            type1_train["type"]="type1"
            type1_test["type"]="type1"
            if debug:    
                print("TYPE 1:")
                print("      - %d molecules in training set."%(len(type1_train)))
                print("      - %d molecules in test set."%(len(type1_test)))
                

            (type2_train,type2_test)=getTrainTestSet("prepared_data/type2.csv",start_date)
            type2_train["type"]="type2"
            type2_test["type"]="type2"

            if debug:  
                print("TYPE 2:")
                print("      - %d molecules in training set."%(len(type2_train)))
                print("      - %d molecules in test set."%(len(type2_test)))
                

            (type1_2_train,type1_2_test)=getTrainTestSet("prepared_data/type1_2.csv",start_date)
            type1_2_train["type"]="type1_2"
            type1_2_test["type"]="type1_2"

            if debug:  
                print("TYPE 1 1/2:")
                print("      - %d molecules in training set."%(len(type1_2_train)))
                print("      - %d molecules in test set."%(len(type1_2_test)))
            
            #smiles_train=list(type1_train.smiles)+list(type2_train.smiles)
            smiles_train=list(type1_train.smiles)+list(type2_train.smiles)+list(type1_2_train.smiles)
            if(len(smiles_train)>0): 
                train=pd.concat([type1_train, type2_train,type1_2_train])
                train_outcome=list(type1_train.type)+list(type2_train.type)+list(type1_2_train.type)
                smiles_test=list(type1_test.smiles)+list(type2_test.smiles)+list(type1_2_test.smiles)
                test=pd.concat([type1_test, type2_test,type1_2_test])
                
                test_outcome=list(type1_test.type)+list(type2_test.type)+list(type1_2_test.type)
                
                morganFpRadius=radius
                fingerprintMethod=AllChem.GetMorganFingerprint


                (fp_train,mol_train,tr_sm)=generateFingerprints(smiles_train,fingerprintMethod,morganFpRadius=morganFpRadius)
                (fp_test,mol_test,te_sm)=generateFingerprints(smiles_test,fingerprintMethod,morganFpRadius=morganFpRadius)
                (predictions,similarity)=getMaxSimilarity(fp_test,fp_train,train_outcome)
                (ranpredictions,ransimilarity)=getRandomPrediction(fp_test,fp_train,train_outcome)

                # print(len(fp_test),len(test["smiles"]),len(fp_test),len(mol_test))
                # print(len(fp_train),len(train["smiles"]),len(fp_train),len(mol_train))
                # print(len(similarity),len(predictions))
                pred=np.array(predictions)
                y_true=np.array(test_outcome)
                mask=pred==y_true
                randomMask=np.array(ranpredictions)==y_true
                #print(mask)
                print("AUC:")
                print(sklearn.metrics.roc_auc_score(mask,similarity))
                for idx,el in enumerate(mask):
                    o.write("{}\t{}\t{}\t{}\t{}\t\"{}\"\t{}\t{}\t{}\n".format(cycle,el,pred[idx],y_true[idx],similarity[idx],smiles_test[idx],randomMask[idx],ransimilarity[idx],start_date))

    o.close()