# conda install scikit-learn
# conda install -c conda-forge scikit-optimize
# conda install -c conda-forge rdkit

from Tools.Clustering.butina import cluster_molecules
from Data.datastructures import Dataset
from Data.utils import read_csv
from Data.Representations.binary_fingerprints import ecfp
import numpy as np

def minlog(x):
    return -np.log10(x)

molecules = read_csv(f"example_data/CHEMBL2047_EC50.csv", smiles_col='smiles', label_col='exp_mean [nM]')

data = Dataset(molecules, transform=ecfp, target_transform=minlog)

data[12]
data.show(12)


cluster_molecules(molecules, cutoff=0.4, scaffold=True)


from Data.Data_prep.splitting import random_split_molecules

train, test, val = random_split_molecules(molecules, test_split=0.2, val_split=0.2)

len(train) / len(molecules)
len(test) / len(molecules)
len(val) / len(molecules)

from MoleculeACE.benchmark import models


from MoleculeACE.benchmark import load_data, models, evaluation, utils

# Setup some variables
dataset = 'CHEMBL287_Ki'
descriptor = utils.Descriptors.ATTENTIVE_GRAPH
algorithm = utils.Algorithms.AFP
path_to_model = 'path_to_model.pkl'



def cross_validate(model, x, y, evaluate, cv: int = 5):
    import numpy as np

    scores = []

    for i in cv:
        mod = model.train(x_train_i, y_train_i)
        pred = mod.predict(x_test_i)
        score = evaluate(pred, y_test_i)
        scores.append(score)

    estimated_score = np.mean(scores)

    return estimated_score

from sklearn.model_selection import StratifiedKFold, ShuffleSplit

# if stratified:
#     has_cliffs = [1 if i in self.cliffs.cliff_mols_soft_consensus else 0 for i in self.smiles_train]
#     skf = StratifiedKFold(n_splits=n_splits, random_state=random_state, shuffle=True)
#     cv_folds = skf.split(self.x_train, has_cliffs)
# else:

ss = ShuffleSplit(n_splits=5, random_state=42)
cv_folds = [(i, j) for i, j in ss.split(data.smiles)]

cv_folds_lst = []
for i in cv_folds:

    pass

data.smiles


def fold_split_random(dataset: Dataset, folds: int = 5, random_state: int = 42):
    from sklearn.model_selection import ShuffleSplit

    ss = ShuffleSplit(n_splits=folds, random_state=random_state)
    return [(i, j) for i, j in ss.split(dataset.smiles)]


def fold_split_knn(dataset, k: int = 10, random_state: int = 42):
    from sklearn.cluster import KMeans

    clust = KMeans(n_clusters=10)
    clust.fit(x)


    pass


arr = []
for idx, m in enumerate(data):

    pass



fold_split_random(data)


[(train_idx, test_idx), (), ()]



## TODO
# make cross-validation class
# cross validation function
# make dataset class
# bayesian optimization class
# evol optimization class


# molecules = read_csv(file)
# train, test = split_mols_butina()     split_mols_knn()
# train_data = Dataset(train)
# test_data = Dataset(test)





if __name__ == '__main__':

    # Load data
    data = load_data(dataset, descriptor=descriptor)

    # Load a model
    model = models.load_model(data, algorithm, path_to_model)
    predictions = model.test_predict()

    # Evaluate your model on activity cliff compounds
    results = evaluation.evaluate(data=data, predictions=predictions)



