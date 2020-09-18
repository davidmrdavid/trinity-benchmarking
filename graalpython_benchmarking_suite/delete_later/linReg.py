import numpy as np

class NormalizedLinearRegression():
    def __init__(self, iterations=20, gamma=0.000001):
        self.gamma = gamma
        self.iterations = iterations

    def fit(self, X, y, w_init=None):
        self.w = w_init if w_init is not None else np.matrix(np.random.randn(X.shape[1], 1))
        for _ in range(self.iterations):
            #a = X * self.w
            #print("step 1")
            #print(type(a))
            #b = a - y
            #print("step 2")
            #tempo = X.T
            #print(".5")
            #print(type(b))
            #print(type(X.T))
            #c = X.T * a # b
            #print("step 3")
            #d = self.gamma * c
            #print("step 4")
            #e = self.w - d
            #print("step 5")
            #self.w -= self.gamma * (X.T * (X * self.w - y)) # TODO: THIS WILL FAIL
            self.w = self.w - self.gamma * (X.T * (X * self.w - y))
        return self

    def predict(self, X):
        try:
            getattr(self, "w")
        except AttributeError:
            raise RuntimeError("You must train the regressor before predicting data!")

        return X * self.w
