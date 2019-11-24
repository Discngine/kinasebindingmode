import pandas as pd
import random
import sys
from rdkit.Chem import AllChem,rdReducedGraphs
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
import multiprocessing
import logging



mtran=np.random.RandomState(seed=3333)

def calc_ergfp( fp1, fp2 ):
    #from https://iwatobipen.wordpress.com/2016/01/16/ergfingerprint-in-rdkit/
    denominator = np.sum( np.dot(fp1,fp1) ) + np.sum( np.dot(fp2,fp2) ) - np.sum( np.dot(fp1,fp2 ))
    numerator = np.sum( np.dot(fp1,fp2) )
    return numerator / denominator

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

def generateFingerprints(smiles,fingerprintMethod,morganFpRadius=2):
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
        
    
    smiles_train=list(type1_train.smiles)+list(type2_train.smiles)
    smiles_test=list(type1_test.smiles)+list(type2_test.smiles)

    n_train_type1=len(type1_train)
    n_train_type2=len(type2_train)
    n_train_type1_2=len(type1_2_train)
    
    (fp_train,mol_train)=generateFingerprints(smiles_train,fingerprintMethod,morganFpRadius=morganFpRadius)
    (fp_test,mol_test)=generateFingerprints(smiles_test,fingerprintMethod,morganFpRadius=morganFpRadius)
    
    resList=[]
    nFound=0
    nCorrectTotal=0
    nNotFound=0
    confusionDict={"type1_ok":0,"type1_false":0,"type2_ok":0,"type2_false":0}

    #now test every molecule in the test set
    for idx,molSmile in enumerate(smiles_test):
        if similarityMethod.__name__=="GetFraggleSimilarity":
            sys.stderr.write(str(idx/len(smiles_test)*100.0)+"\n")
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
                            #print(idx,len(mol_test),len(mol_train),fpidx)
                            
                            similarity,match=similarityMethod(mol_test[idx],mol_train[fpidx])
                            #print(time.time()-t1)
                            #print(similarity)
                        except Exception:
                            #print("Failed for smiles : "+smiles_train[fpidx]+" and "+molSmile)
                            similarity=0.0
                            pass
                    elif similarityMethod.__name__=="TverskySimilarity":
                        similarity=similarityMethod(trainFp,fp1,0.9,0.1)
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

                    if evaluationMode=="probability":
                        
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
                    
                    elif evaluationMode=="top":
                        topSimilarIdx=np.argsort(matchingSimilarities)[-1]
                        predictedBindingMode=matchingBindingModes[topSimilarIdx]
                    
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


def getEvaluationStats(similarityMethod,similarityThreshold,n=2,morganFpRadius=2,evaluationMode="top"):
    results=[]
    for i in range(n):   
       
        if(similarityMethod==calc_ergfp):
            results.append(evaluate(rdReducedGraphs.GetErGFingerprint,similarityMethod,similarityThreshold=similarityThreshold,evaluationMode=evaluationMode))
        else:
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
    #print((similarityMethod.__name__,evaluationMode,morganFpRadius,similarityThreshold,np.mean(ba),np.std(ba),np.mean(f1),np.std(f1),np.mean(mcc),np.std(mcc),np.mean(found),np.std(found),np.mean(notFound),np.std(notFound)))
    print("%s\tevaluationMode:\t%s\tFPRadius\t%d\tsimilarityThreshold:\t%.2f\tba:\t%.5f\t%.5f\tf1:\t%.5f\t%.5f\tmcc\t%.5f\t%.5f\tfound:\t%d\t%.2f\tnotFound\t%d\t%.2f"%(similarityMethod.__name__,evaluationMode,morganFpRadius,similarityThreshold,np.mean(ba),np.std(ba),np.mean(f1),np.std(f1),np.mean(mcc),np.std(mcc),np.mean(found),np.std(found),np.mean(notFound),np.std(notFound)))

if __name__ == "__main__":
    similarityMethods=[DataStructs.TanimotoSimilarity,DataStructs.DiceSimilarity,DataStructs.TverskySimilarity,calc_ergfp]#,GetFraggleSimilarity]
    similarityMethods=[GetFraggleSimilarity]
    fingerprintMethods=[rdReducedGraphs.GetErGFingerprint,AllChem.GetMorganFingerprint]
    threads=[]
    evaluationModes=["top","probability"]
    #similarityMethods=[DataStructs.TanimotoSimilarity,DataStructs.FingerprintSimilarity,DataStructs.DiceSimilarity,DataStructs.AsymmetricSimilarity,DataStructs.BraunBlanquetSimilarity,DataStructs.CosineSimilarity,DataStructs.KulczynskiSimilarity,DataStructs.McConnaugheySimilarity, DataStructs.RogotGoldbergSimilarity, DataStructs.RusselSimilarity,DataStructs.SokalSimilarity,DataStructs.TverskySimilarity]
    
    similarityThresholds=np.arange(0.9,0.25,-0.05)
    for similarityMethod in similarityMethods:
        for similarityThreshold in similarityThresholds:
            for evaluationMode in evaluationModes:
                if similarityMethod != calc_ergfp and similarityMethod != GetFraggleSimilarity:
                    for morganFpRadius in range(1,4):
                        x = multiprocessing.Process(target=getEvaluationStats, args=(similarityMethod,similarityThreshold,10,morganFpRadius,evaluationMode))
                        threads.append(x)
                        x.start()
                        while (len(threads)>=12):
                            for tix,thread in enumerate(threads):
                                if not thread.is_alive():
                                    threads.pop(tix)
                                    break
                            time.sleep(2)
                else : 
                    x = multiprocessing.Process(target=getEvaluationStats, args=(similarityMethod,similarityThreshold,10,0,evaluationMode))
                    threads.append(x)
                    x.start()

                    while (len(threads)>=12):
                        for tix,thread in enumerate(threads):
                            if not thread.is_alive():
                                threads.pop(tix)
                                break
                        time.sleep(2)
        #for index, thread in enumerate(threads):
        #    thread.join()
        
            
        #for acceptanceRatio in acceptanceRatios:
        
                #try:
                #except Exception as err:
                #    print(similarityMethod.__name__+" no results")
                #    print(err)
                #    pass
