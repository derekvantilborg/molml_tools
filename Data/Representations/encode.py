
from Data.Representations import binary_fingerprints
from Data.Representations import classical_descriptors
from Data import Molecule
from rdkit import Chem, DataStructs
import numpy as np


def molecules_to_array(list_of_molecules: list, descriptor: str = 'ecfp'):

    fps = []
    if descriptor == 'ecfp':
        fps = [binary_fingerprints.ecfp(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'maccs':
        fps = [binary_fingerprints.maccs(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'rdkit':
        fps = [binary_fingerprints.rdkit(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'atom_pair':
        fps = [binary_fingerprints.atom_pair(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'topological_torsion':
        fps = [binary_fingerprints.topological_torsion(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'smarts':
        fps = [binary_fingerprints.smarts(m.rdkit) for m in list_of_molecules]
        fps = rdkit_numpy_convert(fps)
    if descriptor == 'drug_like':
        fps = [classical_descriptors.drug_like_descriptor(m.rdkit) for m in list_of_molecules]
        fps = np.array(fps)
    if descriptor == 'whim':
        fps = [classical_descriptors.whim_descriptor(m.rdkit) for m in list_of_molecules]
        fps = np.array(fps)

    return fps


def rdkit_numpy_convert(fp):
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)