# conda install scikit-learn
# conda install -c conda-forge scikit-optimize
# conda install -c conda-forge rdkit
import pandas as pd

# from Tools.Clustering.butina import cluster_molecules
from Data.datastructures import Dataset
from Data.utils import read_csv
from Representations.descriptors import ecfp, maccs
from Representations.strings import smiles_one_hot
from sklearn.ensemble import GradientBoostingRegressor
from Tools.optimize import BayesianOpt
from Tools.metrics import rmse
import numpy as np

def minlog(x):
    return -np.log10(x)

molecules = read_csv(f"example_data/CHEMBL2047_EC50.csv", smiles_col='smiles', label_col='exp_mean [nM]')

data = Dataset(molecules[:50], name='CHEMBL2047', transform=smiles_one_hot, target_transform=minlog)
data.process()

data.show(10)

from Tools.cluster import spectral
from Viz.multivariate import TSNE, PCA
import seaborn as sns

clusters = spectral(molecules, k=10)


tsne = TSNE(n_components=2, perplexity=50, n_iter=500)
tsne.fit(molecules, use_n_principal_components=50)
tsne.show(color_by=clusters, palette=sns.color_palette("hls", 10))

pca = PCA(n_components=2)
pca.fit(molecules)
pca.show(color_by=clusters, palette=sns.color_palette("hls", 10))



from Tools.splitting import stratified_split_molecules

train, test, val = stratified_split_molecules(molecules, labels=clusters)


data = Dataset(molecules, name='CHEMBL2047', transform=ecfp, target_transform=minlog)
data.process()

data.show(13)

hpm = {"learning_rate": [0.1, 0.01],
       "max_depth": [1, 2, 3, 4, 5, 6, 7, 8],
       "n_estimators": [5, 10, 20, 100, 200, 300]}

model = GradientBoostingRegressor

opt = BayesianOpt(model, data)
opt.opt(hpm, rmse, cv=5, n_calls=20)
opt.show()


# def fold_split_knn(dataset, k: int = 10, random_state: int = 42):
#     from sklearn.cluster import KMeans
#
#     clust = KMeans(n_clusters=10)
#     clust.fit(x)


history = [(1,0.7201,0.7201),(2,0.6329,0.6329),(3,0.6305,0.6305),(4,0.6323,0.6305),(5,0.7195,0.6305),(6,0.6137,0.6137),
           (7,0.6201,0.6137),(8,0.6239,0.6137),(9,0.6404,0.6137),(10,0.6264,0.6137),(11,0.6718,0.6137),(12,0.6368,0.6137),
           (13,0.6337,0.6137),(14,0.6502,0.6137),(15,0.6235,0.6137),(16,0.6303,0.6137),(17,0.6171,0.6137),(18,0.6268,0.6137),
           (19,0.6117,0.6117),(20,0.6170,0.6117)]


history = pd.DataFrame( columns=['Iteration', 'Score', 'Best Score'])

history['Score'].tolist()[-1]
len(history['Score'])
pd.DataFrame({'Iteration': [21], 'Score': [0.544], 'Best Score': [0.544]})


## TODO active learning
# split data train test -> make TSNE
# optimize model on train
# train model
# predict on test
# find most uncertain compounds
#



