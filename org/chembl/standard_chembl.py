#!/usr/bin/env python
# coding: utf-8

import rdkit
from rdkit import Chem
import pandas as pd
from molvs import standardize_smiles
import math
import numpy as np
import re
import os
from rdkit.Chem import Descriptors

moe = pd.read_csv('CHEMBL26-chembl_26_activity.csv',  engine='python')
moe = moe[['Molecule ChEMBL ID','Smiles','Standard Value','Standard Units', 'Assay Description']]
moe.head()

old_list = []
smi_list = []
chid_list = []
logs_list = []
dtype_list = []

smi_org = []
logs_org = []

for idx, item in moe.iterrows():
    smi = item['Smiles']
    chid = item['Molecule ChEMBL ID']
    old = smi
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

    mol    = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol

    print (smiles)
    #mol = Chem.AddHs(mol)
    if mol == None: 
        print(smi, 'MolFromSmiles Error')
        continue

    if mol.GetNumAtoms() < 3:
        print(smi, '<100, filterd')
        continue

    m = Chem.MolFromSmiles(smiles)
    weight = Descriptors.MolWt(m)

    value = float(item['Standard Value']);
    print (value)
    unit  = item['Standard Units']
    desc  = item['Assay Description']
   
    if math.isnan(value) : 
        continue
 
    print ("value")
    print (value)
    if math.isclose(value,0) :
        continue

    if unit=='ug.mL-1' or unit=='Ug/ml' or unit=='ug/ml' or unit=='ugA/ml':
        logS = math.log(value/1000.0 / weight, 10)
    elif unit=='mg kg-1':
        logS = math.log(value/1000.0 /weight, 10)
    elif unit=='nM':
        logS = math.log(value/1000.0/1000.0/1000.0, 10) 
    else: 
        continue 
 
    print (logS)
    smi_list.append(smiles)
    chid_list.append(chid)
    logs_list.append(logS)
    dtype_list.append(desc) 
    old_list.append(old)

df = pd.DataFrame({'smiles':smi_list, 'ChEMBL ID':chid_list, 'logS':logs_list, 'description':dtype_list, 'oldsmi':old_list})
df.sort_values(by='smiles')

df.to_csv('chembl_stand.csv', index=False)
