import csv
from rdkit  import Chem  
from data.mol2features import get_available_features_generators, get_features_generator
from data.utils import load_features, save_features
from tqdm import tqdm
import numpy as np

csvPath="stand.csv"
features_generator="rdkit_2d_normalized"
savePath="hello.npz"

print(Chem.__file__)
with open(csvPath) as f:
    reader = csv.reader(f)
    next(reader)  # skip header

    mols = []
    targets = []
    lines = []
    for line in reader:
        smiles = line[0]
        lines.append(line)
        mols.append(Chem.MolFromSmiles(line[0]))
        targets.append(line[1])
    
    print(get_available_features_generators())
    features_generator = get_features_generator(features_generator)
    features_map = map(features_generator, mols)
   
    features = []
    for i, feats in tqdm(enumerate(features_map), total=len(lines)):
    #    print(i)
        features.append(feats)

    print (type(features)) 
    print (type(features[0]))
    save_features(savePath, features)
#    np.savetxt('hello.txt',features)
