import rdkit
import pandas as pd
from molvs import standardize_smiles
import numpy as np
import re
import os
from tqdm.autonotebook import tqdm
import csv
from rdkit import Chem
from rdkit.Chem import PandasTools

ptsdf = PandasTools.LoadSDF('./selected_records.sdf')


sdf = ptsdf[['SMILES','Water solubility Kinetic {measured, converted}', 'pH','Temperature']]
smi_list = []
logs_list = []
pH_list = []
temp_list = []

smi_org = []
logs_org = []

for idx, row in tqdm(sdf.iterrows()):
    smi = row['SMILES']
    value = row['Water solubility Kinetic {measured, converted}'] 
    smi_org.append(smi)
    logs_org.append(value)

#    print (idx, smi)
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
#    print ("standardlized")

    mol = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True) if type(mol) != str else mol

#    print (smiles)
    #mol = Chem.AddHs(mol)
    if mol == None:
        print(smi, 'MolFromSmiles Error')
        continue

    if mol.GetNumAtoms() < 3:
        print(smi, '<100, filterd')
        continue


    smi_list.append(smiles)
    value = row['Water solubility Kinetic {measured, converted}'] 
    logs_list.append(value)
    pH_list.append(row['pH'])
    temp_list.append(row['Temperature']) 

df = pd.DataFrame({'smiles':smi_org, 'logS':logs_org})
df.to_csv('kinect_org.csv', index=False)

df = pd.DataFrame({'smiles':smi_list, 'logS':logs_list, 'pH':pH_list, 'Temperature':temp_list})
df.sort_values(by='smiles')
df.to_csv('kinect_stand.csv', index=False)



