import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from Data.Molecule import Molecule


def ecfp(molecule: Molecule, radius: int = 2, nbits: int = 1024):
    # Calculate the morgan fingerprint
    return AllChem.GetMorganFingerprintAsBitVect(molecule.rdkit, radius, nBits=nbits)


def maccs(molecule: Molecule):
    from rdkit.Chem import MACCSkeys
    return MACCSkeys.GenMACCSKeys(molecule.rdkit)


def rdkit(molecule: Molecule):
    return Chem.RDKFingerprint(molecule.rdkit)


def atom_pair(molecule: Molecule):
    return AllChem.GetHashedAtomPairFingerprintAsBitVect(molecule.rdkit)


def topological_torsion(molecule: Molecule):
    return AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(molecule.rdkit)


def smarts(molecule: Molecule):
    return Chem.PatternFingerprint(molecule.rdkit)

