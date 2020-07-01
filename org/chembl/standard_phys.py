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
moe = moe[['Smiles','Standard Value','Standard Units', 'Assay Description']]
moe.head()

old_list = []
smi_list = []
logs_list = []
dtype_list = []

for idx, item in moe.iterrows():
    smi = item['Smiles']
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

    import pybel
    belmol=pybel.readstring('smi', smiles)
    belmol.make3D() # make 3d coords
    belmol.removeh() 
    xyz=belmol.write('xyz')
    str3d=xyz.split('\n')
    
    print ("xyz")
    positions = []
    for i, atom in enumerate(mol.GetAtoms()):
        tmp = re.split('\s+', str3d[i+2])
        coord = [float(s) for s in tmp[1:]]
        positions.append(coord)

    positions = np.array(positions, dtype = float)  
    array_sum = np.sum(positions)
    array_has_nan = np.isnan(array_sum)
    if array_has_nan == True:
        print(smi, '3D xyz contains NaN')
        continue

    print ("3D")
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
    logs_list.append(logS)
    dtype_list.append(desc) 
    old_list.append(old)

df = pd.DataFrame({'smiles':smi_list, 'logS':logs_list, 'description':dtype_list, 'oldsmi':old_list})
df.sort_values(by='smiles')

df.to_csv('chembl_stand.csv', index=False)
