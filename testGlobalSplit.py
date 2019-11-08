import pandas as pd
import random
import requests as requests
from requests import Session
from requests.auth import HTTPBasicAuth
import sys
import urllib

type1=pd.read_csv("prepared_data/type1.csv")
type1["type"]="type1"
type2=pd.read_csv("prepared_data/type1.csv")
type2["type"]="type2"
type1_2=pd.read_csv("prepared_data/type1.csv")
type1_2["type"]="type1_2"
smiles=list(type1.smiles)+list(type2.smiles)+list(type1_2.smiles)

random.shuffle(smiles)
split=int(len(smiles)/2)
train_data=smiles[:split]
test_data=smiles[split:]
secret=""
host="https://3decision.discngine.cloud/api/v2"
s=Session()
response=s.get(host+'/token',headers = {"x-api-secret":secret})
if response.status_code!=200:
    sys.exit("Could not login")
token=response.json()["data"]["token"]


for test_smiles in test_data[:3]:
    print(test_smiles)
    query=host+"/ligand/search/similarity/smiles/"+urllib.parse.quote(test_smiles,safe='')+"?similarityThreshold=0.8&addConformations=true"
    print(query)
    simSearchResponse=s.get(query,headers={"Authorization":"Bearer "+token})
    print(simSearchResponse)
    if(simSearchResponse.status_code==200):
        print(simSearchResponse.content)