from Data.Molecule import Molecule
from Tools.Clustering.butina import cluster_molecules

import pandas as pd

df = pd.read_csv(f"example_data/CHEMBL2047_EC50.csv")

molecules = [Molecule(smi) for smi in df['smiles']]


cluster_molecules(molecules, cutoff=0.4, scaffold=True)


from Data.Data_prep.splitting import random_split_molecules

train, test, val = random_split_molecules(molecules, test_split=0.2, val_split=0.2)

len(train) / len(molecules)
len(test) / len(molecules)
len(val) / len(molecules)




