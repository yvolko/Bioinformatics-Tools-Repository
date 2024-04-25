from sklearn.base import BaseEstimator


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
            self, n_estimators=10, max_depth=None, max_features=None,
            random_state=SEED
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, x, y):
        self.classes_ = sorted(np.unique(y))

        for i in range(self.n_estimators - 1):
            np.random.seed(self.random_state + i)

            feat_ids = np.random.choice(x.shape[1], self.max_features,
                                        replace=False)
            self.feat_ids_by_tree.append(feat_ids)

            sample_indices = np.random.choice(x.shape[0], x.shape[0],
                                              replace=True)
            x_sampled = x[sample_indices, [feat_ids]]
            y_sampled = y[sample_indices]

            tree = DecisionTreeClassifier(max_depth=self.max_depth,
                                          random_state=self.random_state)
            tree.fit(x_sampled, y_sampled)
            self.trees.append(tree)

        return self

    def predict_proba(self, x):
        probabilities = []

        for i in range(len(self.trees)):
            tree = self.trees[i]
            feat_id = self.feat_ids_by_tree[i]
            x_selected = x[:, feat_id]
            probabilities.append(tree.predict_proba(x_selected))

        return sum(probabilities) / len(probabilities)

    def predict(self, x):
        probas = self.predict_proba(x)
        predictions = np.argmax(probas, axis=1)

        return predictions
