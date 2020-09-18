import numpy as np

class GaussianNMF():
    def __init__(self, iterations=20, components=5):
        self.iterations = iterations
        self.components = components

    def fit(self, X, w_init=None, h_init=None):
        print(w_init)
        print(h_init)
        self.w = w_init if w_init is not None else np.mat(np.random.rand(X.shape[0], self.components))
        self.h = h_init if h_init is not None else np.mat(np.random.rand(self.components, X.shape[1]))
        print("forward")

        """Factorize non-negative matrix."""
        for _ in range(self.iterations):
            print("HERE")
            #print(type(self.w))
            #print(self.w.T)
            print(".....")
            #print(type(X))
            #print(self.w.T * X)
            print("ayayayay")
            a = (self.w.T * X) 
            print(">>>>>>")
            b =  (self.w.T * self.w * self.h) 
            
            print("!!!")
            self.h = np.multiply(self.h, (self.w.T * X) / (self.w.T * self.w * self.h))
            self.w = np.multiply(self.w, (X * self.h.T) / (self.w * (self.h * self.h.T)))
            print("NO...way")

        return self
