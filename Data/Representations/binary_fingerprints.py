import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from Data.Molecule import Molecule


def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)


def ecfp(smiles: str, to_array: bool = True, radius: int = 2, nbits: int = 1024):
    # Calculate the morgan fingerprint
    fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), radius, nBits=nbits)
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp


def maccs(smiles: str, to_array: bool = True):
    from rdkit.Chem import MACCSkeys
    fp = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp

def rdkit(smiles: str, to_array: bool = True):
    fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp

def atom_pair(smiles: str, to_array: bool = True):
    fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(Chem.MolFromSmiles(smiles))
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp

def topological_torsion(smiles: str, to_array: bool = True):
    fp = AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect(Chem.MolFromSmiles(smiles))
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp

def smarts(smiles: str, to_array: bool = True):
    fp = Chem.PatternFingerprint(Chem.MolFromSmiles(smiles))
    if to_array:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return fp
