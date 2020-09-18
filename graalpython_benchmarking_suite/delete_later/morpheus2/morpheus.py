"""
Python package to support 'morpheus' rewrites of linear algebra operations
"""
import copy
import polyglot
import numpy as np
import numpy.core.numeric as N
from numpy.core.numeric import isscalar
from numpy.matrixlib.defmatrix import asmatrix, matrix

class MorpheuInvalidArgsError(TypeError):
    pass

class TensorFromNumpy(object):
    
    def __init__(self,  npArr):
        # TODO: validate!
        self.npArr = npArr

    def unwrap(self):
        return self.npArr

    def rightMatrixMultiplication(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.matmul(self.npArr, otherArr))

    def leftMatrixMultiplication(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.matmul(otherArr, self.npArr))

    def scalarAddition(self, scalar):
        return TensorFromNumpy(self.npArr + scalar)

    def scalarMultiplication(self, scalar):
        return TensorFromNumpy(self.npArr * scalar)

    def crossProduct(self):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.cross(self.npArr, self.npArr.T))

    def rowSum(self):
        return TensorFromNumpy(np.sum(self.npArr, axis=1))

    def columnSum(self):
        return TensorFromNumpy(np.sum(self.npArr, axis=0))

    def elementWiseSum(self):
        return float(np.sum(self.npArr, axis=None))

    def rowWiseAppend(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.vstack([self.npArr, otherArr]))

    def columnWiseAppend(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.hstack([self.npArr, otherArr]))

    def splice(self, begRows, endRows, begCols, endCols):
        return TensorFromNumpy(self.npArr[begRows:endRows+1, begCols:endCols+1]) # TODO: does not use axis, needs better designs

    def matrixAddition(self, otherMatrix):
        if isinstance(otherMatrix, TensorFromNumpy):
            otherMatrix = otherMatrix.npArr
        return TensorFromNumpy(self.npArr + otherMatrix)

    def getLength(self, axis):
        return float(self.npArr.shape[axis])

    def getNumRows(self):
        return self.npArr.shape[0]

    def getNumCols(self):
        return self.npArr.shape[1]

    def scalarExponentiation(self, exp):
        return TensorFromNumpy(np.power(self.npArr, exp))
    
    # TODO: Do we want this?
    def transpose(self):
        return TensorFromNumpy(self.npArr.T)

    # TODO: scalarOps, aggregation, others   

# TODO: this should inherit from np.array, thus behave like one
class NormalizedMatrix(np.ndarray):

    """
    Array functions are created to follow numpy semantics.
    https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html
    
    1. Explicit constructor call
    2. View casting
    3. Creating from new template
    """
    def __array_prepare__(self, obj, context=None):
        pass

    def __array_wrap__(self, out_arr, context=None):
        pass

    def __array_finalize__(self, obj):
        pass

    _SUPPORTED_UFUNCS = {np.add: {1: "__add__", -1: "__radd__"},
                         np.subtract: {1: "__sub__", -1: "__rsub__"},
                         np.divide: {1: "__div__", -1: "__rdiv__"},
                         np.multiply: {1: "__mul__", -1: "__rmul__"},
                         np.power: {1: "__pow__", -1: "__rpow__"}}

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Handle ufunc supported in numpy standard library.
        reference: https://docs.scipy.org/doc/numpy-1.13.0/reference/ufuncs.html

        :param ufunc: ufunc object
        :param method: type of method. In this class, only __call__ is handled
        :param inputs:
        :param kwargs:
        :return: Normalized matrix or matrix or ndarray or numeric
        """
        if ufunc in self._SUPPORTED_UFUNCS and len(inputs) == 2 and method == "__call__":
            order = isinstance(inputs[0], NormalizedMatrix) - isinstance(inputs[1], NormalizedMatrix)
            if order == 1:
                return getattr(inputs[0], self._SUPPORTED_UFUNCS[ufunc][order])(inputs[1], **kwargs)
            if order == -1:
                return getattr(inputs[1], self._SUPPORTED_UFUNCS[ufunc][order])(inputs[0], **kwargs)
            if order == 0 and ufunc is np.multiply:
                return inputs[0].__mul__(inputs[1], **kwargs)
        return NotImplemented


    def __getNormalizedTable(self, tensorS, tensorKs, tensorRs):
        normalizedTable = polyglot.eval(language="thanos",string="")
        normalizedTable.build(tensorS, tensorKs, tensorRs)
        return normalizedTable

    # TODO: 
    # 1. Do need __array_priority__ (?)
    # 2. In MorpheusPy, why does `__getitem__` cause significant performance penalties?
    # 3. Why inherit from matrix and yet call numeric as the super?
    # 
    def __new__(cls, S, Ks, Rs, is_transposed=False):
        obj = N.ndarray.__new__(cls, None)
        #print("Test 1")
        tensorS = TensorFromNumpy(S)
        tensorKs = list(map(lambda x: TensorFromNumpy(x), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x), Rs))
        obj.is_transposed = is_transposed
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs
        #print("In __new__")
        obj.normalizedTable = obj.__getNormalizedTable(tensorS, tensorKs, tensorRs)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        S = getattr(obj, 'S', None)
        tensorS = TensorFromNumpy(S)
        Ks = getattr(obj, 'Ks', None) 
        Rs = getattr(obj, 'Rs', None)
        tensorKs = list(map(lambda x: TensorFromNumpy(x), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x), Rs))
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs
        obj.is_tranposed = getattr(obj, 'is_transposed', None)
        obj.normalizedTable = obj.__getNormalizedTable(tensorS, tensorKs, tensorRs)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        #TODO: provide more information
        return "NormalizedMatrix"

    def __getitem__(self, index):
        # TODO: provide better error
        return NotImplementedError

    """ columnSum, rowSum, elementWiseSum """
    def sum(self, axis=None, dtype=None, out=None):
        if axis == 0:
            return self.normalizedTable.columnSum().unwrap()
        elif axis == 1:
            return self.normalizedTable.rowSum().unwrap()
        return self.normalizedTable.elementWiseSum().unwrap()

    """ scalarExponentiation """
    def __pow__(self, other):
        #print("IN POW")
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
        normMatrix.normalizedTable = self.normalizedTable.scalarExponentiation(other)
        return normMatrix
    def __ipow__(self, other):
        #print("IN POW2")
        output = self.__pow__(other)
        self.normalizedTable = output.normalizedTable #TODO: confirm this is fine
        return self
    def __rpow__(self, other):
        #print("IN POW3")
        return self.__pow__(other)

    """ matrixAddition and scalarAddition"""
    def __add__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarAddition(other)
            return normMatrix
        return self.normalizedTable.matrixAddition(other).unwrap()
    def __iadd__(self, other):
        output =  self.__add__(other) 
        self.normalizedTable = output.normalizedTable
        return self
    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)
    def __isub__(self, other):
        return self.__iadd__(-other)
    def __rsub__(self, other):
        return self.__radd__(-other)

    
    """ scalarMultiplication, crossProd, LMM and RMM """
    #TODO: morpheusPy also handles the case where 'other' has no
    # __rmul__ method, but that's for another day.
    def __mul__(self, other):
        print("NO FUCKING WAY")
        #print("HI")
        ##print(self.is_transposed)
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix
        if isinstance(other, NormalizedMatrix):
            #TODO: stamp1!

            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented
        if isinstance(other, (N.ndarray, list,tuple)):
           #print("BRANCH gOOD")
           #if self.is_transposed:
           #    print("do RMM")
           #    return self.normalizedTable.rightMatrixMultiplication(TensorFromNumpy(other)).unwrap()
           #else:
           return self.normalizedTable.leftMatrixMultiplication(TensorFromNumpy(other)).unwrap()
        return NotImplemented

    def __rmul__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix

        if isinstance(other, NormalizedMatrix):
            #TODO: stamp!

            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented
        if isinstance(other, (N.ndarray, list,tuple)):
           #if self.is_transposed:
           #    return self.normalizedTable.leftMatrixMultiplication(TensorFromNumpy(other)).unwrap()
           #else:
           return self.normalizedTable.rightMatrixMultiplication(TensorFromNumpy(other)).unwrap()
        return NotImplemented

    def __imul__(self, other):
        if not isscalar(other):
            return NotImplemented

        self.normalizedTable = self.normalizedTable.scalarMultiplication(other)
        return self

    def __truediv__(self, other):
        return self.__mul__(1/other)
    def __rtruediv__(self, other):
        return self.__rdiv__(1/other)
    def __itruediv__(self, other):
        return self.__idiv__(1/other)

    @property
    def T(self):
        #print("HERE")
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
        normMatrix.is_transposed = not(self.is_transposed)
        normMatrix.normalizedTable = self.normalizedTable.transpose()
        return normMatrix
        #normMatrix = self.normalizedTable.transpose()
        #newNormalizedTable = copy.deepcopy(self)#NormalizedMatrix()
        #newNormalizedTable.normalizedTable = normMatrix
        #return newNormalizedTable


# TODO: this should inherit from np.array, thus behave like one
class NormalizedMatrix2222(np.ndarray):

    """
    Array functions are created to follow numpy semantics.
    https://docs.scipy.org/doc/numpy-1.13.0/user/basics.subclassing.html
    
    1. Explicit constructor call
    2. View casting
    3. Creating from new template
    """
    def __array_prepare__(self, obj, context=None):
        pass

    def __array_wrap__(self, out_arr, context=None):
        pass

    def __array_finalize__(self, obj):
        pass

    _SUPPORTED_UFUNCS = {np.add: {1: "__add__", -1: "__radd__"},
                         np.subtract: {1: "__sub__", -1: "__rsub__"},
                         np.divide: {1: "__div__", -1: "__rdiv__"},
                         np.multiply: {1: "__mul__", -1: "__rmul__"},
                         np.power: {1: "__pow__", -1: "__rpow__"}}

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Handle ufunc supported in numpy standard library.
        reference: https://docs.scipy.org/doc/numpy-1.13.0/reference/ufuncs.html

        :param ufunc: ufunc object
        :param method: type of method. In this class, only __call__ is handled
        :param inputs:
        :param kwargs:
        :return: Normalized matrix or matrix or ndarray or numeric
        """
        if ufunc in self._SUPPORTED_UFUNCS and len(inputs) == 2 and method == "__call__":
            order = isinstance(inputs[0], NormalizedMatrix) - isinstance(inputs[1], NormalizedMatrix)
            if order == 1:
                return getattr(inputs[0], self._SUPPORTED_UFUNCS[ufunc][order])(inputs[1], **kwargs)
            if order == -1:
                return getattr(inputs[1], self._SUPPORTED_UFUNCS[ufunc][order])(inputs[0], **kwargs)
            if order == 0 and ufunc is np.multiply:
                return inputs[0].__mul__(inputs[1], **kwargs)
        return NotImplemented


    def __getNormalizedTable(self, tensorS, tensorKs, tensorRs):
        normalizedTable = polyglot.eval(language="thanos",string="")
        normalizedTable.build(tensorS, tensorKs, tensorRs)
        return normalizedTable

    # TODO: 
    # 1. Do need __array_priority__ (?)
    # 2. In MorpheusPy, why does `__getitem__` cause significant performance penalties?
    # 3. Why inherit from matrix and yet call numeric as the super?
    # 
    def __new__(cls, S, Ks, Rs, is_transposed=False):
        obj = N.ndarray.__new__(cls, None)
        #print("Test 1")
        tensorS = TensorFromNumpy(S)
        tensorKs = list(map(lambda x: TensorFromNumpy(x), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x), Rs))
        obj.is_transposed = is_transposed
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs
        #print("In __new__")
        obj.normalizedTable = obj.__getNormalizedTable(tensorS, tensorKs, tensorRs)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        S = getattr(obj, 'S', None)
        tensorS = TensorFromNumpy(S)
        Ks = getattr(obj, 'Ks', None) 
        Rs = getattr(obj, 'Rs', None)
        tensorKs = list(map(lambda x: TensorFromNumpy(x), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x), Rs))
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs
        obj.is_tranposed = getattr(obj, 'is_transposed', None)
        obj.normalizedTable = obj.__getNormalizedTable(tensorS, tensorKs, tensorRs)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        #TODO: provide more information
        return "NormalizedMatrix"

    def __getitem__(self, index):
        # TODO: provide better error
        return NotImplementedError

    """ columnSum, rowSum, elementWiseSum """
    def sum(self, axis=None, dtype=None, out=None):
        if axis == 0:
            return self.normalizedTable.columnSum().unwrap()
        elif axis == 1:
            return self.normalizedTable.rowSum().unwrap()
        return self.normalizedTable.elementWiseSum().unwrap()

    """ scalarExponentiation """
    def __pow__(self, other):
        #print("IN POW")
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
        normMatrix.normalizedTable = self.normalizedTable.scalarExponentiation(other)
        return normMatrix
    def __ipow__(self, other):
        #print("IN POW2")
        output = self.__pow__(other)
        self.normalizedTable = output.normalizedTable #TODO: confirm this is fine
        return self
    def __rpow__(self, other):
        #print("IN POW3")
        return self.__pow__(other)

    """ matrixAddition and scalarAddition"""
    def __add__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarAddition(other)
            return normMatrix
        return self.normalizedTable.matrixAddition(other).unwrap()
    def __iadd__(self, other):
        output =  self.__add__(other) 
        self.normalizedTable = output.normalizedTable
        return self
    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-other)
    def __isub__(self, other):
        return self.__iadd__(-other)
    def __rsub__(self, other):
        return self.__radd__(-other)

    
    """ scalarMultiplication, crossProd, LMM and RMM """
    #TODO: morpheusPy also handles the case where 'other' has no
    # __rmul__ method, but that's for another day.
    def __mul__(self, other):
        print("NO FUCKING WAY")
        #print("HI")
        ##print(self.is_transposed)
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix
        if isinstance(other, NormalizedMatrix):
            #TODO: stamp1!

            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented
        if isinstance(other, (N.ndarray, list,tuple)):
           #print("BRANCH gOOD")
           #if self.is_transposed:
           #    print("do RMM")
           #    return self.normalizedTable.rightMatrixMultiplication(TensorFromNumpy(other)).unwrap()
           #else:
           return self.normalizedTable.leftMatrixMultiplication(TensorFromNumpy(other)).unwrap()
        return NotImplemented

    def __rmul__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix

        if isinstance(other, NormalizedMatrix):
            #TODO: stamp!

            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented
        if isinstance(other, (N.ndarray, list,tuple)):
           #if self.is_transposed:
           #    return self.normalizedTable.leftMatrixMultiplication(TensorFromNumpy(other)).unwrap()
           #else:
           return self.normalizedTable.rightMatrixMultiplication(TensorFromNumpy(other)).unwrap()
        return NotImplemented

    def __imul__(self, other):
        if not isscalar(other):
            return NotImplemented

        self.normalizedTable = self.normalizedTable.scalarMultiplication(other)
        return self

    def __truediv__(self, other):
        return self.__mul__(1/other)
    def __rtruediv__(self, other):
        return self.__rdiv__(1/other)
    def __itruediv__(self, other):
        return self.__idiv__(1/other)

    @property
    def T(self):
        #print("HERE")
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
        normMatrix.is_transposed = not(self.is_transposed)
        normMatrix.normalizedTable = self.normalizedTable.transpose()
        return normMatrix
        #normMatrix = self.normalizedTable.transpose()
        #newNormalizedTable = copy.deepcopy(self)#NormalizedMatrix()
        #newNormalizedTable.normalizedTable = normMatrix
        #return newNormalizedTable
