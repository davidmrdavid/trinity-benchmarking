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

class NormalizedLogisticRegression:
    def __init__(self, iterations=20, gamma=0.000001):
        self.gamma = gamma
        self.iterations = iterations

    def fit(self, X, y, w_init=None):
        self.w = w_init if w_init is not None else np.matrix(np.random.randn(X.shape[1], 1))
        for i in range(self.iterations):
            print(type(self.w))
            print(type(X))
            print(X.isMorpheus)
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            bb = X * self.w
            cc = self.w - self.gamma * (y / (1 + 2.71 ** bb))
            raise NotImplementedError
            #self.w = self.w - self.gamma * (X.T * (y / (1 + 2.71 ** (X * self.w))))
        return self

    def predict(self, X):
        try:
            getattr(self, "w")
        except AttributeError:
            raise RuntimeError("You must train the regressor before predicting data!")

        return X * self.w
