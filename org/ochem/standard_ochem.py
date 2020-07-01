#!/usr/bin/env python
# coding: utf-8

#find run the following shell script to preprocess the csv file
#cat ../orig_data/ochem.csv | grep -v "#6\]" | grep -v "#7\]" | sort | uniq > ochem2.csv
#cat ochem/ochem_train.csv | sort | uniq | grep -v "\[#6\]" | grep -v "\[#7\]" > ochem2.csv
#tile smiles,logS up

import rdkit
from rdkit import Chem
import pandas as pd
from molvs import standardize_smiles
from tqdm.autonotebook import tqdm
import numpy as np
import re
import os

moe = pd.read_csv('standtmp.csv', engine='python')
moe = moe[['SMILES','Water solubility {measured, converted}','UNIT {Water solubility}', 'Dataset', 'Temperature', 'UNIT {Temperature}', 'pH','UNIT {pH}', 'Quality code', 'UNIT {Quality code}']]

moe.head()

smi_list = []
logs_list = []
dtype_list = []

smi_org = []
logs_org = []

#for idx, item in tqdm(moe.iterrows()):
for idx, item in moe.iterrows():
    smi = item['SMILES']
    value = str(float(item['Water solubility {measured, converted}']) * -1.0)
    smi_org.append(smi)
    logs_org.append(value)

    print (idx, smi)
    if smi is np.nan:
        print(smi, 'smiles Error')
        continue
    if len(smi) < 2:
        print(smi, 'smiles too short')
        continue

    try:
        smiles = standardize_smiles(smi)
    except rdkit.Chem.rdchem.AtomValenceException:
        print ("Error in try")
        continue;
    except rdkit.Chem.rdchem.AtomKekulizeException:
        print ("Error in try")
        continue;
    except rdkit.Chem.rdchem.KekulizeException:
        print ("Error in try")
        continue;
    else:
        print ("smiles Pass")

    if smi is np.nan:
        print(smi, 'standardize Error')
        continue    
    print ("standardlized")

    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol

    print (smiles)
    #mol = Chem.AddHs(mol)
    if mol == None: 
        print(smi, 'MolFromSmiles Error')
        continue

    if mol.GetNumAtoms() < 3:
        print(smi, '<100, filterd')
        continue

    smi_list.append(smiles)
    value = str(float(item['Water solubility {measured, converted}']) * -1.0)
    print (item['Water solubility {measured, converted}'])
    print (value)
    logs_list.append(value)
    dtype_list.append(item['Dataset']) 


df = pd.DataFrame({'smiles':smi_org, 'logS':logs_org})
df.to_csv('ochem_org.csv', index=False)

df = pd.DataFrame({'smiles':smi_list, 'logS':logs_list, 'dtype':dtype_list})
df.sort_values(by='smiles')

df.to_csv('ochem_stand.csv', index=False)
