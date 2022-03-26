""" Molecule class that holds everything about a molecule. Based on RDkit"""
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from Data.Representations.binary_fingerprints import ecfp

class Molecule:
    def __init__(self, smiles: str, y: float = None, id: str = None):
        self.smiles = smiles
        self.y = y
        self.id = id

    def mol(self):
        return Chem.MolFromSmiles(self.smiles)

    def ecfp(self, to_array: bool = False, radius: int = 2, nbits: int = 1024):
        # Calculate the morgan fingerprint
        return ecfp(self.smiles, to_array, radius, nbits)

    def show(self, size: tuple = (500, 500), kekulize: bool = True):
        """ Plot an image of the molecule """
        from rdkit.Chem import Draw
        Draw.ShowMol(self.mol(), size=size, kekulize=kekulize)

    def __repr__(self):
        return self.smiles


class Dataset:
    def __init__(self, molecules: list = None, smiles: list = None, ids: list = None, ys: list = None,
                 transform=None, target_transform=None):

        self.molecules = molecules

        # if molecules is not None:
        #     self.smiles, self.ids, self.ys = [], [], []
        #     for m in molecules:
        #         self.smiles.append(m.smiles)
        #         self.ids.append(m.id)
        #         self.ys.append(m.y)
        # else:
        #     self.smiles, self.ids, self.ys = smiles, ids, ys

        self.transform = transform
        self.target_transform = target_transform
        self.filename = None

    #
    # class molecules():
    #     def __init__(self, smiles, ids=None, ys=None):
    #         self.smiles = smiles
    #         self.ids = ids
    #         self.ys = ys
    #
    #     def molecule(self, idx):
    #         smi = self.smiles[idx]
    #         id = self.ids[idx] if self.ids is not None else None
    #         y = self.ys[idx] if self.ys is not None else None
    #         return Molecule(smi, y, id)
    #
    #     def __getitem__(self, idx):
    #         super.__init__()
    #         return self.smiles[idx]

    def show(self, idx, size: tuple = (500, 500), kekulize: bool = True):
        if self.molecules is not None:
            self.molecules[idx].show(size, kekulize)

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):

        x = self.molecules[idx].smiles
        y = self.molecules[idx].y
        if self.transform:
            x = self.transform(x)
        if self.target_transform:
            y = self.target_transform(y)
        return x, y




#
# import numpy as np
# #
# def minlog(x):
#     return -np.log10(x)
#
# def maccs(smi):
#     from rdkit.Chem import MACCSkeys
#     return MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smi))
#
# molecules = read_csv(f"example_data/CHEMBL2047_EC50.csv", smiles_col='smiles', label_col='exp_mean [nM]')
#
# data = Dataset(molecules, transform=maccs, target_transform=minlog)
# # data.read_csv(f"example_data/CHEMBL2047_EC50.csv", smiles_col='smiles', label_col='exp_mean [nM]')



# molecules[25].show()
#

# data.molecules[12]

# # # data.show(12)
# data.show(25)
#
# data[25]

