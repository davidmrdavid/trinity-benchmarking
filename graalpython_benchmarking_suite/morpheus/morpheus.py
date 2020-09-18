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
    
    def __init__(self,  npArr, foreign_backend=False):
        self.npArr = np.matrix(npArr)
        self.foreign_backend = foreign_backend

    def unwrap(self):
        if self.foreign_backend:
            return self.npArr.tolist()
        return self.npArr

    def removeAbstractions(self):
        return self.unwrap()

    def rightMatrixMultiplication(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.matmul(self.npArr, otherArr), self.foreign_backend)

    def leftMatrixMultiplication(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.matmul(otherArr, self.npArr), self.foreign_backend)

    def scalarAddition(self, scalar):
        return TensorFromNumpy(self.npArr + scalar, self.foreign_backend)

    def scalarMultiplication(self, scalar):
        return TensorFromNumpy(self.npArr * scalar, self.foreign_backend)

    def crossProduct(self):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.cross(self.npArr, self.npArr.T), self.foreign_backend)

    def rowSum(self):
        return TensorFromNumpy(np.sum(self.npArr, axis=1), self.foreign_backend)

    def columnSum(self):
        return TensorFromNumpy(np.sum(self.npArr, axis=0), self.foreign_backend)

    def elementWiseSum(self):
        return float(np.sum(self.npArr, axis=None))

    def rowWiseAppend(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.vstack([self.npArr, otherArr]), self.foreign_backend)

    def columnWiseAppend(self, otherArr):
        if isinstance(otherArr, TensorFromNumpy):
            otherArr = otherArr.npArr
        return TensorFromNumpy(np.hstack([self.npArr, otherArr]), self.foreign_backend)

    def splice(self, begRows, endRows, begCols, endCols):
        return TensorFromNumpy(self.npArr[begRows:endRows+1, begCols:endCols+1], self.foreign_backend)

    def matrixAddition(self, otherMatrix):
        if isinstance(otherMatrix, TensorFromNumpy):
            otherMatrix = otherMatrix.npArr
        return TensorFromNumpy(self.npArr + otherMatrix, self.foreign_backend)

    def getLength(self, axis):
        return float(self.npArr.shape[axis])

    def getNumRows(self):
        return self.npArr.shape[0]

    def getNumCols(self):
        return self.npArr.shape[1]

    def scalarExponentiation(self, exp):
        return TensorFromNumpy(np.power(self.npArr, exp), self.foreign_backend)
    
    def transpose(self):
        return TensorFromNumpy(self.npArr.T, self.foreign_backend)

# TODO: this should inherit from np.array, thus behave like one
class NormalizedMatrix(matrix): #(np.array): << SWITCH 2 DIS IF THINGS BREAK

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


    def __getNormalizedTable(self, tensorS, tensorKs, tensorRs, Sempty):
        normalizedTable = polyglot.eval(language="morpheusDSL",string="")
        # TODO: Check if S is empty
        normalizedTable.build(tensorS, tensorKs, tensorRs, Sempty)
        return normalizedTable

    # TODO: 
    # 1. Do need __array_priority__ (?)
    # 2. In MorpheusPy, why does `__getitem__` cause significant performance penalties?
    # 3. Why inherit from matrix and yet call numeric as the super?
    # 
    def __new__(cls, S, Ks, Rs, is_transposed=False, foreign_backend=False, Scheckable=True):
        obj = N.ndarray.__new__(cls, None)
        obj.foreign_backend = foreign_backend
        tensorS = TensorFromNumpy(S, obj.foreign_backend)
        tensorKs = list(map(lambda x: TensorFromNumpy(x, obj.foreign_backend), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x, obj.foreign_backend), Rs))

        obj.is_transposed = is_transposed
        obj.Scheckable = Scheckable
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs

        Sempty = False
        if Scheckable and len(S) == 0:
            Sempty = True
        if str(type(S)) == "<class 'foreign'>":
            obj.normalizedTable = obj.__getNormalizedTable(S, Ks, Rs, Sempty)
            obj.foreign_backend = True
        else:
            obj.normalizedTable = obj.__getNormalizedTable(tensorS, tensorKs, tensorRs, Sempty)
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return

        obj.foreign_backend = getattr(obj, 'foreign_backend', None)
        S = getattr(obj, 'S', None)
        tensorS = TensorFromNumpy(S, obj.foreign_backend)
        Ks = getattr(obj, 'Ks', None) 
        Rs = getattr(obj, 'Rs', None)

        tensorKs = list(map(lambda x: TensorFromNumpy(x, obj.foreign_backend), Ks))
        tensorRs = list(map(lambda x: TensorFromNumpy(x, obj.foreign_backend), Rs))
        obj.S = S
        obj.Ks = Ks
        obj.Rs = Rs
        obj.is_tranposed = getattr(obj, 'is_transposed', None)
        obj.Scheckable = getattr(obj, 'Scheckable', None)
        if str(type(S)) == "<class 'foreign'>":
            obj.normalizedTable = obj.__getNormalizedTable(S, Ks, Rs)
            obj.foreign_backend = True
        else:
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
        output = None
        if axis == 0:
            output = self.normalizedTable.columnSum()
        elif axis == 1:
            output = self.normalizedTable.rowSum()
        else:
            output = self.normalizedTable.elementWiseSum()
        if self.foreign_backend:
            # TODO: DO BETTER, THIS IS A  HACK KNOWING BACK IS R
            return np.matrix(output.unwrap()).reshape(output.getNumRows()[0], output.getNumCols()[0])
        return output.unwrap()

    """ scalarExponentiation """
    def __pow__(self, other):
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
        normMatrix.normalizedTable = self.normalizedTable.scalarExponentiation(other)
        return normMatrix
    def __ipow__(self, other):
        output = self.__pow__(other)
        self.normalizedTable = output.normalizedTable # TODO: confirm this is fine
        return self
    def __rpow__(self, other):
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
        if str(type(S)) == "<class 'foreign'>": # if True:
            output = self.matrixAddition(other.scalarMultiplication(-1))
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs, Scheckable=self.Scheckable)
            normMatrix.normalizedTable = output
            return normMatrix
        return self.__add__(-other)
    def __isub__(self, other):
        return self.__iadd__(-other)
    def __rsub__(self, other):
        return self.__radd__(-other)

    
    """ scalarMultiplication, crossProd, LMM and RMM """
    #TODO: morpheusPy also handles the case where 'other' has no
    # __rmul__ method, but that's for another day.
    def __mul__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix
        if isinstance(other, NormalizedMatrix):
            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented

        if isinstance(other, Faux):
            output = self.normalizedTable.leftMatrixMultiplication(other.obj)
            faux = Faux(output)
            return faux
        if isinstance(other, (N.ndarray, list,tuple)):
           output = self.normalizedTable.leftMatrixMultiplication(TensorFromNumpy(other, self.foreign_backend))
           if self.foreign_backend:
               got = list(output.unwrap())
               ncol = output.getNumCols()[0] # TODO: fix
               nrow = output.getNumRows()[0] # TODO: fix
               return np.matrix(got).reshape(nrow, ncol)
           else:
               return output.unwrap()

        elif str(type(self.S)) == "<class 'foreign'>": # if True:
            output = self.normalizedTable.leftMatrixMultiplication(other)
            faux = Faux(output)
            return faux
        return NotImplemented

    def __rmul__(self, other):
        if isscalar(other):
            normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs)
            normMatrix.normalizedTable = self.normalizedTable.scalarMultiplication(other)
            return normMatrix

        if isinstance(other, NormalizedMatrix):
            if self.stamp == other.stamp and self.is_transposed ^ other.is_transposed:
                return self.normalizedTable.crossProduct().unwrap()
            else:
                # TODO: can we implement this?
                return NotImplemented
        if isinstance(other, (N.ndarray, list,tuple)):
           output = self.normalizedTable.rightMatrixMultiplication(TensorFromNumpy(other))
           if self.foreign_backend:
               got = list(output.unwrap())   # TODO: why list() ?
               ncol = output.getNumCols()[0] # TODO: fix
               nrow = output.getNumRows()[0] # TODO: fix
               return np.matrix(got).reshape(nrow, ncol)
           else:
               return output.unwrap()
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
        normMatrix = NormalizedMatrix(self.S, self.Ks, self.Rs, foreign_backend=self.foreign_backend, Scheckable=self.Scheckable)
        normMatrix.is_transposed = not(self.is_transposed)
        normMatrix.normalizedTable = self.normalizedTable.transpose()
        return normMatrix

# ==================

class Faux(matrix): #(np.array): << SWITCH 2 DIS IF THINGS BREAK

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

    # TODO: 
    # 1. Do need __array_priority__ (?)
    # 2. In MorpheusPy, why does `__getitem__` cause significant performance penalties?
    # 3. Why inherit from matrix and yet call numeric as the super?
    # 
    def __new__(cls, mat):
        obj = N.ndarray.__new__(cls, None)
        obj.obj = mat
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        mm = getattr(obj, 'obj', None)
        obj.obj = mm

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
        output = None
        if axis == 0:
            output = self.obj.columnSum()
        elif axis == 1:
            output = self.obj.rowSum()
        else:
            output = self.obj.elementWiseSum()
            return output
        if True:
            return np.matrix(output.unwrap()).reshape(output.getNumRows()[0], output.getNumCols()[0]) # TODO: DO BETTER, THIS IS A  HACK KNOWING BACK IS R
        return output.unwrap()

    """ scalarExponentiation """
    def __pow__(self, other, expo):
        mm = Faux(other.obj)
        mm.obj = other.obj.scalarExponentiation(self)
        return mm
    def __ipow__(self, other, expo):
        output = self.__pow__(other, expo)
        other.obj = output.obj #TODO: confirm this is fine
        return other
    def __rpow__(self, other, expo):
        return self.__pow__(other, expo)

    """ matrixAddition and scalarAddition"""
    def __add__(self, other):
        if isscalar(other):
            mm = Faux(self.obj)
            mm.obj = self.obj.scalarAddition(other)
            return mm
        return self.obj.matrixAddition(other).unwrap()

    def __radd__(self, other):
        return self.__add__(other)
    def __iadd__(self, other):
        output =  self.__add__(other) 
        self.obj = output.obj
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Faux):
            faux = Faux(self.obj.matrixAddition(other.obj.scalarMultiplication(-1)))
            return faux
        faux = Faux(self.obj.matrixAddition(other.scalarMultiplication(-1)))
        return faux
    def __isub__(self, other):
        return self.__iadd__(-other)
    def __rsub__(self, other):
        faux = Faux(self.obj.matrixAddition(other.scalarMultiplication(-1)))
        return faux

    
    """ scalarMultiplication, crossProd, LMM and RMM """
    #TODO: LMM AND RMM ARE STIL SWTICHEEEEED
    #TODO: morpheusPy also handles the case where 'other' has no
    # __rmul__ method, but that's for another day.
    def __mul__(self, other):
        if isscalar(other):
            mm = Faux(self.obj)
            mm.obj = self.obj.scalarMultiplication(other)
            return mm
        if True: #isinstance(other, (N.ndarray, list,tuple)): #TODO
           if not(isinstance(other, (N.ndarray, list, tuple))):
               output = self.obj.rightMatrixMultiplication(other)
               mm = Faux(output)
               return mm
           else:
               output = self.obj.rightMatrixMultiplication(other.obj)
               mm = Faux(output)
               return mm
        return NotImplemented

    def __rmul__(self, other):
        if isscalar(other):
            mm = Faux(self.obj)
            mm.obj = self.obj.scalarMultiplication(other)
            return mm

        if True: #isinstance(other, (N.ndarray, list,tuple)): # TODO
           if not(isinstance(other, (N.ndarray, list, tuple))):
               output = self.obj.leftMatrixMultiplication(other)
               mm = Faux(output)
               return mm
           else:
               output = self.obj.leftMatrixMultiplication(other.obj)
               mm = Faux(output)
               return mm
        return NotImplemented

    def __imul__(self, other):
        if not isscalar(other):
            return NotImplemented

        self.obj = self.obj.scalarMultiplication(other)
        return self

    def __truediv__(self, other):
        #x = other.invert()
        faux = Faux(self.obj)
        faux.obj = faux.obj.divisionArr(other.obj)
        return faux
    def __rtruediv__(self, other):
        faux = Faux(self.obj)
        faux.obj = faux.obj.divisionArr(other)
        faux.obj = faux.obj.invert()
        return faux
    def __itruediv__(self, other):
        return NotImplemented

    @property
    def T(self):
        mm = Faux(self.obj)
        mm.obj = mm.obj.transpose()
        return mm
