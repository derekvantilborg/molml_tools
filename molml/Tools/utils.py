
from molml.Datastructures.dataset import Dataset
from molml.Tools.splitting import fold_split_random
from typing import Callable
import keras as K


def cross_validate(model: Callable, dataset: Dataset, evaluate: Callable, cv: int = 5, random_state: int = 42,
                   verbose: bool = True, to_array: bool = True):
    """ Cross validate a model in folds

    Args:
        model: (Callable) Any model that has a .fit(x, y) method and a .predict(x) method
        dataset: (Dataset) A Dataset object that holds molecules
        evaluate: (Callable) Any function that evaluates and scores (true, pred). Should work with numpy arrays
        cv: (int/list[tuple]) Takes an iterator of tuples of train/test indices. Create random folds if int is given
        random_state: (int) random seed when no cv is given
        verbose: (bool) print progress yes/no

    Returns: Estimated cross-validation score

    """
    import numpy as np
    from copy import deepcopy

    if type(cv) is int:
        cv = fold_split_random(dataset, random_state=random_state, folds=cv)

    x_train = dataset.get_x(to_array=to_array)
    y_train = dataset.get_y(to_array=to_array)
    x_test = dataset.get_x(to_array=to_array)
    y_test = dataset.get_y(to_array=to_array)
    mod = None

    scores = []
    for fold, indices in enumerate(cv):
        # Copy a version of the model for this cv-fold

        mod = deepcopy(model)

        # Get fold data
        train_idx, test_idx = indices
        x_train_ = np.array([x_train[i] for i in train_idx])
        y_train_ = np.array([y_train[i] for i in train_idx])
        x_test_ = np.array([x_test[i] for i in test_idx])
        y_test_ = np.array([y_test[i] for i in test_idx])

        # Train a model
        mod.fit(x_train_, y_train_)

        # predict test data
        pred = mod.predict(x_test_)

        # Evaluate
        score = evaluate(y_test_, pred)
        scores.append(score)
        if verbose:
            print(f"fold {fold+1}/{len(cv)}: {score}")

        K.clear_session()
        del mod

    K.clear_session()
    del model

    estimated_score = np.mean(scores)

    if verbose:
        print(f"estimated score: {estimated_score}")

    return estimated_score


def create_search_space(hyperparameters: dict):
    """ Convert a dict of lists into a hyperparameter search space for skopt"""
    from skopt.space import Categorical
    search_space = []
    for k, v in hyperparameters.items():
        search_space.append(Categorical(categories=v, name=k))

    return search_space
