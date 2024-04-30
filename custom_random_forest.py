import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from concurrent.futures import ThreadPoolExecutor


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
            self, random_state=21, n_estimators=10,
            max_depth=None, max_features=None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, x, y, n_jobs):
        self.classes_ = sorted(np.unique(y))

        def train_tree(seed):
            np.random.seed(seed)

            feat_ids = np.random.choice(x.shape[1], self.max_features,
                                        replace=False)
            self.feat_ids_by_tree.append(feat_ids)

            sample_indices = np.random.choice(x.shape[0], x.shape[0],
                                              replace=True)

            x_sampled = x[sample_indices, feat_ids]
            y_sampled = y[sample_indices]

            if x_sampled.ndim == 1:
                x_sampled = x_sampled.reshape(-1, 1)

            tree = DecisionTreeClassifier(max_depth=self.max_depth,
                                          random_state=self.random_state)
            tree.fit(x_sampled, y_sampled)
            self.trees.append(tree)

        seeds = [self.random_state + i for i in range(self.n_estimators - 1)]

        with ThreadPoolExecutor(n_jobs) as pool:
            pool.map(train_tree, seeds)

        return self

    def predict_proba(self, x, n_jobs):
        probabilities = []

        def compute_tree_proba(tree, feat_id, x_selected):
            if x_selected.ndim == 1:
                x_selected = x_selected.reshape(-1, 1)
            return tree.predict_proba(x_selected)

        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            results = [pool.submit(compute_tree_proba, self.trees[i],
                                   self.feat_ids_by_tree[i],
                                   x[:, feat_id]) for i, feat_id in
                       enumerate(self.feat_ids_by_tree)]

        for result in results:
            probabilities.append(result.result())

        return np.mean(probabilities, axis=0)

    def predict(self, x, n_jobs):
        probas = self.predict_proba(x, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
