import numpy as np
from Data.Molecule import Molecule


def drug_like_descriptor(molecule: Molecule):
    from rdkit.Chem.Descriptors import ExactMolWt
    from rdkit.Chem import Descriptors, rdMolDescriptors, rdmolops, QED, Crippen, rdchem

    # https://sharifsuliman1.medium.com/understanding-drug-likeness-filters-with-rdkit-and-exploring-the-withdrawn-database-ebd6b8b2921e

    weight = ExactMolWt(molecule.rdkit)
    logp = Descriptors.MolLogP(molecule.rdkit)
    h_bond_donor = Descriptors.NumHDonors(molecule.rdkit)
    h_bond_acceptors = Descriptors.NumHAcceptors(molecule.rdkit)
    rotatable_bonds = Descriptors.NumRotatableBonds(molecule.rdkit)
    atoms = rdchem.Mol.GetNumAtoms(molecule.rdkit)
    heavy_atoms = rdchem.Mol.GetNumHeavyAtoms(molecule.rdkit)
    molar_refractivity = Crippen.MolMR(molecule.rdkit)
    topological_polar_surface_area = QED.properties(molecule.rdkit).PSA
    formal_charge = rdmolops.GetFormalCharge(molecule.rdkit)
    rings = rdMolDescriptors.CalcNumRings(molecule.rdkit)

    descr = np.array([weight, logp, h_bond_donor, h_bond_acceptors, rotatable_bonds, atoms, heavy_atoms,
                      molar_refractivity, topological_polar_surface_area, formal_charge, rings])
    return descr


def whim_descriptor(molecule: Molecule, seed=0xf00d):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

    # Add hydrogens
    mh = Chem.AddHs(molecule.rdkit)
    # compute conformers with the experimental torsion knowledge distance geometry method
    # https://pubs.acs.org/doi/10.1021/acs.jcim.5b00522
    ps = AllChem.ETKDGv2()
    ps.randomSeed = seed
    # Use distance geometry to obtain initial coordinates for a molecule
    AllChem.EmbedMolecule(mh, ps)
    # calculate WHIM 3D descriptor
    whim = rdMolDescriptors.CalcWHIM(mh)

    return np.array(whim)
