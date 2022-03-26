import numpy as np
from Data.Molecule import Molecule
from rdkit import Chem

def drug_like_descriptor(smiles: str):
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem import Descriptors, rdMolDescriptors, rdmolops, QED, Crippen, rdchem

    # https://sharifsuliman1.medium.com/understanding-drug-likeness-filters-with-rdkit-and-exploring-the-withdrawn-database-ebd6b8b2921e
    mol = Chem.MolFromSmiles(smiles)
    weight = ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_bond_donor = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    atoms = rdchem.Mol.GetNumAtoms(mol)
    heavy_atoms = rdchem.Mol.GetNumHeavyAtoms(mol)
    molar_refractivity = Crippen.MolMR(mol)
    topological_polar_surface_area = QED.properties(mol).PSA
    formal_charge = rdmolops.GetFormalCharge(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)

    descr = np.array([weight, logp, h_bond_donor, h_bond_acceptors, rotatable_bonds, atoms, heavy_atoms,
                      molar_refractivity, topological_polar_surface_area, formal_charge, rings])
    return descr


def whim_descriptor(molecule: Molecule, seed=0xf00d):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    # Add hydrogens
    mh = Chem.AddHs(mol)
    # compute conformers with the experimental torsion knowledge distance geometry method
    # https://pubs.acs.org/doi/10.1021/acs.jcim.5b00522
    ps = AllChem.ETKDGv2()
    ps.randomSeed = seed
    # Use distance geometry to obtain initial coordinates for a molecule
    AllChem.EmbedMolecule(mh, ps)
    # calculate WHIM 3D descriptor
    whim = rdMolDescriptors.CalcWHIM(mh)

    return np.array(whim)
