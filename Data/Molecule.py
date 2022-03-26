""" Molecule class that holds everything about a molecule. Based on RDkit"""
from rdkit import Chem


class Molecule:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.rdkit = Chem.MolFromSmiles(smiles)
        self.scaffold = None

    def show(self, size: tuple = (500, 500), kekulize: bool = True):
        """ Plot an image of the molecule """
        from rdkit.Chem import Draw
        Draw.ShowMol(self.rdkit, size=size, kekulize=kekulize)

    def __repr__(self):
        return self.smiles
