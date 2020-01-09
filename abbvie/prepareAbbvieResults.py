import pandas as pd
import numpy as np 
from rdkit.Chem import AllChem
from rdkit import Chem
from sklearn.metrics import balanced_accuracy_score, f1_score, matthews_corrcoef
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('max_colwidth', -1)
debug=0
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')    
import time

from rdkit import DataStructs
#print(abv_type1_fp)

def predictExternalSet(input_fp):
    response=[]
    for ix,abv_fp in enumerate(input_fp):
        maxSim=0.0
        curRefIx=-1
        if(abv_fp!=''):
            for ref_ix,ref_fp in enumerate(fp_train):
                if(ref_fp!=''):
                    sim=DataStructs.TverskySimilarity(abv_fp,ref_fp,0.9,0.1)
                    if(sim>maxSim):
                        curRefIx=ref_ix
                        maxSim=sim
            predictedBindingMode=getBindingModeFromSmiles(smiles_train[curRefIx])
            if(maxSim>0.3):
                response.append(predictedBindingMode)
            else :
                response.append(0)
        else: 
            response.append(0)
    return response

def getBindingModeFromSmiles(smiles):
        if(sum(type1_train.smiles==smiles) or sum(type1_test.smiles==smiles)): return "type1"
        if(sum(type2_train.smiles==smiles) or sum(type2_test.smiles==smiles)): return "type2"
        if(sum(type1_2_train.smiles==smiles) or sum(type1_2_test.smiles==smiles)): return "type1_2"


def getTrainTestSet(fileName,fraction,random_seed=None):
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


data=pd.read_csv("/Users/peter/Desktop/AbbvieKinConformations.csv",sep=";")
print(len(np.unique(data["pdb"])),"\tunique pdb structures")
print(np.sum(data["errorType"]!="none"),"\twith some sort of error")
print("DFG conformations:")
print(data["DFG_Type"].value_counts())
print("alpha C helix conformations:")
print(data["Helix_Type"].value_counts())

abv_type1=data.iloc[np.where((~data["smiles"].isnull()) & (data["smiles"]!=" ") & (data["DFG_Type"]=="in") & (data["Helix_Type"]=="in") )[0],:]
abv_type1_2=data.iloc[np.where((~data["smiles"].isnull()) & (data["DFG_Type"]=="in") & (data["Helix_Type"]=="out"))[0],:]
abv_type2=data.iloc[np.where((~data["smiles"].isnull()) & (data["smiles"]!="") & (data["DFG_Type"]=="out") & (data["Helix_Type"]=="out") & (data["back"]==True))[0],:]

print("Type I:", len(abv_type1))
print("Type I1/2:", len(abv_type1_2))
#print("Type I1/2 out-like:", len(np.where((~data["smiles"].isnull()) & (data["smiles"]!="") & (data["DFG_Type"]=="in") & (data["Helix_Type"]=="out-like"))[0]))
print("Type II:", len(abv_type2))
#print("Type II (ligand or not):", len(np.where((data["DFG_Type"]=="out") & (data["Helix_Type"]=="out"))[0]))
#print((data["smiles"]).isnull())
print(abv_type2["smiles"].to_string())



for cycle in range(10):

    (type1_train,type1_test)=getTrainTestSet("../prepared_data/type1.csv",0.5)
    type1_train["type"]="type1"
    type1_test["type"]="type1"
    if debug:    
        print("TYPE 1:")
        print("      - %d molecules in training set."%(len(type1_train)))
        print("      - %d molecules in test set."%(len(type1_test)))
        

    (type2_train,type2_test)=getTrainTestSet("../prepared_data/type2.csv",0.5)
    type2_train["type"]="type2"
    type2_test["type"]="type2"

    if debug:  
        print("TYPE 2:")
        print("      - %d molecules in training set."%(len(type2_train)))
        print("      - %d molecules in test set."%(len(type2_test)))
        

    (type1_2_train,type1_2_test)=getTrainTestSet("../prepared_data/type1_2.csv",0.5)
    type1_2_train["type"]="type1_2"
    type1_2_test["type"]="type1_2"

    if debug:  
        print("TYPE 1 1/2:")
        print("      - %d molecules in training set."%(len(type1_2_train)))
        print("      - %d molecules in test set."%(len(type1_2_test)))
        

    smiles_train=list(type1_train.smiles)+list(type2_train.smiles)
    smiles_train=list(type1_train.smiles)+list(type1_test.smiles)+list(type2_test.smiles)+list(type2_train.smiles)
    smiles_test=list(type1_test.smiles)+list(type2_test.smiles)

    n_train_type1=len(type1_train)
    n_train_type2=len(type2_train)
    n_train_type1_2=len(type1_2_train)

    morganFpRadius=3
    fingerprintMethod=AllChem.GetMorganFingerprint


    (fp_train,mol_train)=generateFingerprints(smiles_train,fingerprintMethod,morganFpRadius=morganFpRadius)
    (abv_type1_fp,abv_type1_mol)=generateFingerprints(abv_type1["smiles"],fingerprintMethod,morganFpRadius=morganFpRadius)
    (abv_type1_2_fp,abv_type1_2_mol)=generateFingerprints(abv_type1_2["smiles"],fingerprintMethod,morganFpRadius=morganFpRadius)
    (abv_type2_fp,abv_type2_mol)=generateFingerprints(abv_type2["smiles"],fingerprintMethod,morganFpRadius=morganFpRadius)


    responseType1=np.array(predictExternalSet(abv_type1_fp))

    notFoundType1=np.sum(responseType1=='0')
    #print(notFoundType1," Not found type 1")
    responseType1=responseType1[responseType1!='0']
    #print(np.isnan(np.array(responseType1)))
    respType1=responseType1=="type1"
    #print(np.array(responseType1)=="type1")
    #print(np.sum(np.array(responseType1)=="type1")," TP")
    #print(np.sum(np.array(responseType1)=="type2")," FN")

    responseType2=np.array(predictExternalSet(abv_type2_fp))

    notFoundType2=np.sum(responseType2=='0')
    #print(notFoundType2," Not found type 2")

    responseType2=responseType2[responseType2!='0']
    respType2=responseType2=="type2"
    #print(np.array(responseType2)=="type2")
    #print(np.sum(np.array(responseType2)=="type2")," TN")
    #print(np.sum(np.array(responseType2)=="type1")," FP")

    y_true=np.hstack((np.repeat(1,len(respType1)),np.repeat(0,len(respType2))))
    y_predict=np.hstack((respType1,respType2))

    #print(y_true)
    #print(y_predict)
    ba=balanced_accuracy_score(y_true,y_predict)
    f1=f1_score(y_true,y_predict)
    mcc=matthews_corrcoef(y_true,y_predict)

    print(ba,f1,mcc)