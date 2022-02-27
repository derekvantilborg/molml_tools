from Data.Molecule import Molecule


def clusterfp(fps, cutoff=0.4):
    from rdkit.ML.Cluster import Butina
    from rdkit import DataStructs
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def cluster_molecules(molecules, cutoff=0.4, scaffold=True):
    """ Cluster smiles based on their Murcko scaffold using the Butina algorithm:

    D Butina 'Unsupervised Database Clustering Based on Daylight's Fingerprint and Tanimoto Similarity:
    A Fast and Automated Way to Cluster Small and Large Data Sets', JCICS, 39, 747-750 (1999)
    """
    from Data.Representations import binary_fingerprints

    if scaffold:
        from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
        from rdkit.Chem import AllChem
        # Make Murcko scaffolds
        scaffolds = [GetScaffoldForMol(mol.rdkit) for mol in molecules]
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in scaffolds]
    else:
        fps = [binary_fingerprints.ecfp(mol) for mol in molecules]

    # Cluster fingerprints
    clusters = clusterfp(fps, cutoff=cutoff)
    # smi_clust = [tuple([smiles[idx] for idx in c]) for c in clusters]

    mol_cluster = []
    for clust in clusters:
        mol_cluster.append([molecules[idx] for idx in clust])

    return mol_cluster


