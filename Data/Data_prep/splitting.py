import random
random_seed = 42


def random_split_molecules(molecules, test_split: float = 0.2, val_split: float = 0.1):
    """"""

    if sum([test_split, val_split]) > 1:
        print('Sum of test + val splits cannot be higher than 1. This leaves no training data')

    train, test, val = [], [], []
    random.seed(random_seed)

    for mol in molecules:
        random_nr = random.uniform(0, 1)
        if random_nr < val_split:
            val.append(mol)
        elif val_split < random_nr < (test_split + val_split):
            test.append(mol)
        else:
            train.append(mol)

    return train, test, val


def stratified_split_molecules():
    pass


def cluster_split_molecules():
    pass


def cluster_stratified_split_molecules():
    pass

def fold_split_random():
    pass

def fold_split_stratified():
    pass

def fold_split_knn():
    pass

def fold_split_butina():
    pass

