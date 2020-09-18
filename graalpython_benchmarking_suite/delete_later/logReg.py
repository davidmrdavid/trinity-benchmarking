# Copyright 2018 Side Li and Arun Kumar
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
#from sklearn.base import BaseEstimator, RegressorMixin
#TODO: is inplace substract not workjing??

class NormalizedLogisticRegression: #BaseEstimator, RegressorMixin):
    def __init__(self, iterations=20, gamma=0.000001):
        self.gamma = gamma
        self.iterations = iterations

    def fit(self, X, y, w_init=None):
        self.w = w_init if w_init is not None else np.matrix(np.random.randn(X.shape[1], 1))
        #print("S1")
        #s1 = X * self.w
        #print("S2")
        #s2 = 2.71 ** s1
        #print("S3")
        #print(type(s2))
        #s3 = 1 + s2
        #print("S4")
        #print(type(y))
        #print(type(s3))
        #s4 = y / s3
        #print("S5")
        #print(type(s4))
        #print(type(s4.obj))
        #print(type(X.T))
        #s5 = X.T * s4
        #print("S6")
        #s6 = self.gamma * s5
        #print("GOOOO")
        #exit(-1) 
        #print(self.iterations) 
        for i in range(self.iterations):
            #print(i)
            #self.w = self.w - self.gamma * (X.T * (y / (1 + np.exp(X * self.w))))
            self.w = self.w - self.gamma * (X.T * (y / (1 + 2.71 ** (X * self.w))))
        return self

    def predict(self, X):
        try:
            getattr(self, "w")
        except AttributeError:
            raise RuntimeError("You must train the regressor before predicting data!")

        return X * self.w
