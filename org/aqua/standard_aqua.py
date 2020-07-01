#!/usr/bin/env python
# coding: utf-8

from rdkit import Chem
import pandas as pd
from molvs import standardize_smiles
from tqdm.autonotebook import tqdm
import numpy as np
import re
import os

moe = pd.read_csv('aqua.csv')

moe.head()

smi_list = []
logs_list = []
dtype_list = []

smi_org = []
logs_org = []

for idx, item in moe.iterrows():
    smi   = item.smiles
    value = item.logS
    smi_org.append(smi)
    logs_org.append(value)
    
    print(idx, smi)
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
    value = item.logS
    logs_list.append(value)

df = pd.DataFrame({'smiles':smi_org, 'logS':logs_org})
df.to_csv('aqua_org.csv', index=False)


df = pd.DataFrame({'smiles':smi_list, 'logS':logs_list})
df.sort_values(by='smiles')
df.to_csv('aqua_stand.csv', index=False)