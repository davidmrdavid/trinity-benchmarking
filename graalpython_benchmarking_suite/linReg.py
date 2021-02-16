import numpy as np
from time import time

class NormalizedLinearRegression():
    def __init__(self, iterations=20, gamma=0.000001):
        self.gamma = gamma
        self.iterations = iterations

    def fit(self, X, y, w_init=None):
        self.w = w_init if w_init is not None else np.matrix(np.random.randn(X.shape[1], 1))
        
        for _ in range(self.iterations):
            #bb = (X * self.w)
            #print("......", type(bb), bb.dimm())
            #bb = bb - y
            #print("_______", type(bb), bb.dimm())
            #cc = X.T
            #print("=======", type(cc))
            #dd = (cc * bb)
            #print("^^^^^^^", type(dd), dd.dimm())
            #ee = self.w - dd
            
            #raise UU
            #start = int(time() * 1000)
            tmp = X * self.w
            #end = int(time() * 1000)
            #print("Step 1:", str(end - start))
            
            #start = int(time() * 1000)
            tmp = tmp - y
            #end = int(time() * 1000)
            #print("Step 2:", str(end - start))
            
#             start = int(time() * 1000)
#             XT = X.T
#             end = int(time() * 1000)
#             print("T:", str(end - start))
            
            
            #start = int(time() * 1000)
            tmp = (tmp.T * X).T
            #end = int(time() * 1000)
            #print("Step 3:", str(end - start))
            
            #start = int(time() * 1000)
            tmp = self.gamma * tmp
            #end = int(time() * 1000)
            #print("Step 4:", str(end - start))
            
            #start = int(time() * 1000)
            self.w = self.w - tmp
            #end = int(time() * 1000)
            #print("Step 5:", str(end - start))
            #self.w = self.w - self.gamma * (X.T * (X * self.w - y))
        return self

    def predict(self, X):
        try:
            getattr(self, "w")
        except AttributeError:
            raise RuntimeError("You must train the regressor before predicting data!")

        return X * self.w
