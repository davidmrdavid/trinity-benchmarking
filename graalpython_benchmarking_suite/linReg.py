import numpy as np

class NormalizedLinearRegression():
    def __init__(self, iterations=20, gamma=0.000001):
        self.gamma = gamma
        self.iterations = iterations

    def fit(self, X, y, w_init=None):
        self.w = w_init if w_init is not None else np.matrix(np.random.randn(X.shape[1], 1))
        for _ in range(self.iterations):
            self.w = self.w - self.gamma * (X.T * (X * self.w - y))
        return self

    def predict(self, X):
        try:
            getattr(self, "w")
        except AttributeError:
            raise RuntimeError("You must train the regressor before predicting data!")

        return X * self.w
