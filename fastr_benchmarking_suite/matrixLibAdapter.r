library(Matrix)

MatrixLibAdapter2 <- setClass(
  "MatrixLibAdapter2",
   slot = c(
     rowSum = "ANY",
     columnSum = "ANY",
     leftMatrixMultiplication = "ANY",
     rightMatrixMultiplication = "ANY",
     columnWiseAppend = "ANY",
     transpose = "ANY",
     crossProduct = "ANY",
     crossProductDuo = "ANY",
     elementWiseSqrt = "ANY",
     rowWiseAppend = "ANY",
     diagonal = "ANY",
     splice = "ANY",
     matrixAddition = "ANY",
     getNumRows = "ANY",
     getNumCols = "ANY",
     scalarAddition = "ANY",
     scalarExponentiation = "ANY",
     elementWiseSum = "ANY"
   )
)


setMethod("initialize", "MatrixLibAdapter2",

    function(.Object) {
	.Object@rowSum = function(x) {
            z <- rowSum(x)
            return(z)
        }
        .Object@columnSum = function(x) {
            z <- columnSum(x)
            return(z)
        }
        .Object@leftMatrixMultiplication = function(x, y) {
            z <- leftMatrixMultiplication(x, y)
            if(TRUE){ #SPARSE check
                z <- as.matrix(z)
            }
            return(z)
        }
        .Object@scalarAddition = function(x, y){
            return(scalarAddition(x, y))
        }
        .Object@scalarExponentiation = function(x, y){
            return(scalarExponentiation(x, y))
        }
        .Object@rightMatrixMultiplication = function(x, y) {
            z <- rightMatrixMultiplication(x, y)
            if(TRUE){ #SPARSE check
                z <- as.matrix(z)
            }
            return(z)
        }
        .Object@columnWiseAppend = function(x, y) {
            z <- columnWiseAppend(x, y)
            return(z)
        }
        .Object@transpose = function(x) {
            z <- transpose(x)
            return(z)
        }
        .Object@crossProduct = function(x) {
            z <- crossProduct(x)
            return(z)
        }
        .Object@crossProductDuo = function(x, y) {
            z <- crossProductDuo(x, y)
            return(z)
        }
        .Object@elementWiseSum = function(x) {
            z <- elementWiseSum(x)
            return(z)
        }
	.Object@elementWiseSqrt = function(x) {
            z <- elementWiseSqrt(x)
            return(z)
        }
        .Object@rowWiseAppend = function(x, y) {
            z <- rowWiseAppend(x, y)
            return(z)
        }
        .Object@diagonal = function(x) {
            z <- diagonal(x)
            return(z)
        }
        .Object@splice = function(x, rowBeg, rowEnd, colBeg, colEnd) {
            z <- splice(x, rowBeg, rowEnd, colBeg, colEnd)
            return(z);
        }

        .Object@matrixAddition = function(x, y) {
            z <- matrixAddition(x, y)
	    return(z);

        }
        .Object@getNumCols = function(x) {
            z <- getNumCols(x)
            return(z);


        }
        .Object@getNumRows = function(x) {
            z <- getNumRows(x)
            return(z);

        }
        return(.Object)
})



#setMethod("initialize", "MatrixLibAdapter2",
#  function(.Object) {
#    .Object@columnSum = function(x){
#      return(columnSum(x))
#    }
#    return(.Object)
#  }
#)

MatrixLibAdapter <- setClass(
  "MatrixLibAdapter",
  slot = c(
    matrix = "ANY",
    foreignBackend = "ANY",
    scalarAddition = "ANY",
    scalarExponentiation = "ANY",
    rightMatrixMultiplication = "ANY",
    leftMatrixMultiplication = "ANY",
    scalarMultiplication = "ANY",
    crossProduct = "ANY",
    crossProductDuo = "ANY",
    rowSum = "ANY",
    columnSum = "ANY",
    elementWiseSum = "ANY",
    elementWiseSqrt = "ANY",
    rowWiseAppend = "ANY",
    columnWiseAppend = "ANY",
    matrixAddition = "ANY",
    transpose = "ANY",
    splice = "ANY",
    getNumRows = "ANY",
    getNumCols = "ANY",
    unwrap = "ANY",
    removeAbstractions = "ANY",
    invert = "ANY",
    divisionArr = "ANY",
    diagonal = "ANY"
  )
)

setMethod("initialize", "MatrixLibAdapter",

    function(.Object, matrix, foreignBackend=FALSE) {
        print("BAD NEWS")

        #TODO: This is quite boilerplate-y. Can we do better?
        .Object@matrix <- matrix
        .Object@foreignBackend <- foreignBackend
        .Object@scalarAddition = function(x) {
            return(scalarAddition(.Object, x))
        }
        .Object@scalarExponentiation = function(x){
            return(scalarExponentiation(.Object, x))
        }
        .Object@rightMatrixMultiplication = function(x) {
            return(rightMatrixMultiplication(.Object, x))
        }
        .Object@leftMatrixMultiplication = function(x) {
            return(leftMatrixMultiplication(.Object, x))
        }
        .Object@scalarMultiplication = function(x) {
            return(scalarMultiplication(.Object, x))
        }
        .Object@crossProduct = function() {
            return(crossProduct(.Object, x))
        }
        .Object@crossProductDuo = function(x) {
            return(crossProductDuo(.Object, x))
        }
        .Object@rowSum = function() {
            return(rowSum(.Object))
        }
        .Object@columnSum = function() {
            return(columnSum(.Object))
        }
        .Object@elementWiseSum = function() {
            return(elementWiseSum(.Object))
        }
        .Object@elementWiseSqrt = function() {
            return(elementWiseSqrt(.Object));
        }
        .Object@rowWiseAppend = function(x) {
            return(rowWiseAppend(.Object, x)) 
        }
        .Object@columnWiseAppend = function(x) {
            return(columnWiseAppend(.Object, x))
        }
        .Object@matrixAddition = function(x) {
            return(matrixAddition(.Object,x ))
        }
        .Object@transpose = function() {
            return(transpose(.Object))
        }
        .Object@splice = function(rowBeg, rowEnd, colBeg, colEnd) {
            return(splice(.Object, rowBeg, rowEnd, colBeg, colEnd))
        }
        .Object@getNumCols = function() {
            return(getNumCols(.Object))
        }
        .Object@getNumRows = function() {
            return(getNumRows(.Object))
        }
        .Object@removeAbstractions = function() {
            return(removeAbstractions(.Object))
        }
        .Object@unwrap = function() { # make sure unwrap is no longer used
            return(removeAbstractions(.Object))
        }
        .Object@invert = function() {
            return(invert(.Object))
        }
        .Object@divisionArr = function(x) {
            return(divisionArr(.Object, x))
        }
        .Object@diagonal = function() {
            return(diagonal(.Object))
        }
        return(.Object)
})


setMethod("show", "MatrixLibAdapter", function(object){
    print(object@matrix)
})

# MatrixLibAdapter methods ====================================================
# TODO: how to handle the @matrix redirection boilerplate?


# RMM
# TODO: this could be simplified further
setGeneric("rightMatrixMultiplication", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- tensor %*% otherMatrix
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "ANY", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt = FALSE) {

        result <- rightMatrixMultiplication(tensor@matrix, otherMatrix, FALSE); #foreignBackendOpt = tensor@foreignBackend)
        return(result)
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "polyglot.value", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt = FALSE) {

	#otherMatrixStrict = matrix(otherMatrix$removeAbstractions(), otherMatrix$getNumRows(), otherMatrix$getNumCols())
        otherMatrixStrict = otherMatrix$removeAbstractions()
        result <- rightMatrixMultiplication(tensor@matrix, otherMatrixStrict, FALSE);#foreignBackendOpt = TRUE)
        return(result)
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "numeric", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt = FALSE) {
        result <- rightMatrixMultiplication(tensor@matrix, as.matrix(otherMatrix), FALSE);#foreignBackendOpt = tensor@foreignBackend)
        return(result)
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "matrix", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt = FALSE) {
        result <- rightMatrixMultiplication(tensor@matrix, otherMatrix, FALSE);#foreignBackendOpt = tensor@foreignBackend)
        return(result);
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "dgCMatrix", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt = FALSE) {
        result <- rightMatrixMultiplication(tensor@matrix, otherMatrix, FALSE);#foreignBackendOpt = tensor@foriegnBackend)
        return(result);
})

setMethod("rightMatrixMultiplication", c("MatrixLibAdapter", "MatrixLibAdapter", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- rightMatrixMultiplication(tensor@matrix, otherMatrix@matrix, FALSE );#, foreignBackendOpt=tensor@foreignBackend || otherMatrix@foreignBackend)
        return(result);
})

# 2. LMM

setGeneric("leftMatrixMultiplication", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {

    result <- otherMatrix %*% tensor
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor)
    return(result)
})

setMethod("leftMatrixMultiplication", c("MatrixLibAdapter", "matrix", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- leftMatrixMultiplication(tensor@matrix, otherMatrix, foreignBackendOpt)
        #result@foreignBackend = tensor@foreignBackend
        return(result);
})

setMethod("leftMatrixMultiplication", c("MatrixLibAdapter", "polyglot.value", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
	#otherMatrixStrict = matrix(otherMatrix$removeAbstractions(), otherMatrix$getNumRows(), otherMatrix$getNumCols())
        otherMatrixStrict = otherMatrix$removeAbstractions()
        result <- leftMatrixMultiplication(tensor@matrix, otherMatrixStrict, foreignBackendOpt)
	#result@foreignBackend = TRUE
        return(result);
})

setMethod("leftMatrixMultiplication", c("MatrixLibAdapter", "dgCMatrix", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- leftMatrixMultiplication(tensor@matrix, otherMatrix, foreignBackendOpt)
        #result@foreignBackend = tensor@foreignBackend
        return(result);
})

setMethod("leftMatrixMultiplication", c("MatrixLibAdapter", "MatrixLibAdapter", "ANY"), 
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        #foreignBackendOpt <- tensor@foreignBackend || otherMatrix@foreignBackend
        result <- leftMatrixMultiplication(tensor@matrix, otherMatrix@matrix, foreignBackendOpt)
        return(result);
})

# Scalar Addition

setGeneric("scalarAddition", function(tensor, otherMatrix) {
    result <- tensor + otherMatrix
    #newTensor <- MatrixLibAdapter(matrix=result)
    #return(newTensor);
    return(result);
})

setMethod("scalarAddition", c("MatrixLibAdapter", "numeric"), 
    function(tensor, otherMatrix) {
        result <- scalarAddition(tensor@matrix, otherMatrix)
        #result@foreignBackend = tensor@foreignBackend
        return(result);
})

setGeneric("scalarExponentiation", function(tensor, number) {
        result <- tensor ^ number #number ^ tensor, would fail (?)
        newTensor <- MatrixLibAdapter(matrix=result)
        return(newTensor)
})

setMethod("scalarExponentiation", c("MatrixLibAdapter", "numeric"),
    function(tensor, number) {
        result <- scalarExponentiation(tensor@matrix, number)
        #result@foreignBackend = tensor@foreignBackend
        return(result)
})

# Scalar Multiplication

setGeneric("scalarMultiplication", function(tensor, scalar) {
    result <- tensor * scalar
    #newTensor <- MatrixLibAdapter(matrix=result)
    #return(newTensor);
    return(result);
})

setMethod("scalarMultiplication", c("MatrixLibAdapter", "numeric"),
    function(tensor, scalar) {
        result <- scalarMultiplication(tensor@matrix, scalar)
        #result@foreignBackend = tensor@foreignBackend
        return(result);
})

# Cross Product

setGeneric("crossProduct", function(tensor) {
    result <- tensor * normMatrix
    #newTensor <- MatrixLibAdapter(matrix=result)
    #return(newTensor);
    return(result);
})

setMethod("crossProduct", c("MatrixLibAdapter"), function(tensor) {
    result <- crossProduct(tensor@matrix)
    #result@foreignBackend = tensor@foreignBackend
    return(result);
})

# Row Sum

setGeneric("rowSum", function(tensor, foreignBackendOpt=FALSE) {

    result <- rowSums(tensor)
    newTensor <- MatrixLibAdapter(matrix=result, foreignBackendOpt)
    return(newTensor);
})

setMethod("rowSum", c("MatrixLibAdapter", "ANY"), function(tensor, foreignBackendOpt=FALSE) {
    result <- rowSum(tensor@matrix, foreignBackendOpt=tensor@foreignBackend)
    return(result);
})


# Column Sum

setGeneric("columnSum", function(tensor, foreignBackendOpt=FALSE) {
    #tStart <- as.numeric(Sys.time())*1000;
    result <- colSums(tensor)
    #tEnd <- as.numeric(Sys.time())*1000;
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #total <- tEnd - tStart
    #print(sprintf("COLSUM: %d | %d x %d", total, nrow(tensor), ncol(tensor)));

    return(result);
})

setMethod("columnSum", c("MatrixLibAdapter"), function(tensor, foreignBackendOpt=FALSE) {
    result <- columnSum(tensor@matrix, FALSE)#foreignBackendOpt=tensor@foreignBackend)
    return(result);
})

# Element-wise Sum

setGeneric("elementWiseSum", function(tensor) {
    result <- sum(tensor)
    return(result)
})

setMethod("elementWiseSum", c("MatrixLibAdapter"), function(tensor) {
    result <- elementWiseSum(tensor@matrix)
    return(result);
})

# Row-wise Append

setGeneric("rowWiseAppend", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- rbind(tensor, otherMatrix)
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("rowWiseAppend", c("MatrixLibAdapter", "MatrixLibAdapter", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- rowWiseAppend(tensor@matrix, otherMatrix@matrix, FALSE);#foreignBackendOpt=tensor@foreignBackend)
        return(result);
})

setMethod("rowWiseAppend", c("MatrixLibAdapter", "matrix", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- rowWiseAppend(tensor@matrix, otherMatrix, FALSE);#foreignBackendOpt=tensor@foreignBackend)
        return(result);
})

# Column-wise Append
setGeneric("columnWiseAppend", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    # TODO: this is a little hacky, we need a more robust way of checking types
    result <- NULL;
    #if(is.null(dim(tensor)) && is.null(dim(otherMatrix))) {
    #    print("HERE")
    #    result <- append(tensor, otherMatrix)
    #}
    #else {
    #    print("OR HERE")
        result <- cbind(tensor, otherMatrix)
    #}
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("columnWiseAppend", c("MatrixLibAdapter", "ANY", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- columnWiseAppend(tensor@matrix, otherMatrix, FALSE);#foreignBackendOpt)
        #result@foreignBackend <- tensor@foreignBackend
        return(result);
})

setMethod("columnWiseAppend", c("ANY", "MatrixLibAdapter", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- columnWiseAppend(tensor, otherMatrix@matrix, FALSE);#foreignBackendOpt)
        #result@foreignBackend <- tensor@foreignBackend
        return(result);
})

setMethod("columnWiseAppend", c("MatrixLibAdapter", "MatrixLibAdapter", "ANY"),
    function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
        result <- columnWiseAppend(tensor@matrix, otherMatrix@matrix, FALSE);#)foreignBackendOpt = (tensor@foreignBackend || otherMatrix@foreignBackend))
        return(result);
})

# Matrix Addition

setGeneric("matrixAddition", function(tensor, otherMatrix) {
    result <- tensor + otherMatrix#otherMatrix@matrix
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=otherMatrix@foreignBackend || tensor@foreignBackend)
    #newTensor@foreignBackend <- otherMatrix@foreignBackend || tensor@foreignBackend # OR THIS ONE. IF EITHER, REALLY
    #return(newTensor);
    return(result);
})

# getNumRows

setGeneric("getNumRows", function(tensor) {
    numRows <- nrow(tensor)
    return(numRows)# as.numeric(numRows[1]));
})

setMethod("getNumRows", c("MatrixLibAdapter"), function(tensor) {
    result <- getNumRows(tensor@matrix)
    return(result);
})

# getNumCols

setGeneric("getNumCols", function(tensor) {
    numCols <- ncol(tensor)
    return(numCols);
})

setMethod("getNumCols", c("MatrixLibAdapter"), function(tensor) {
    numCols <- getNumCols(tensor@matrix)
    return(numCols);
})

# transpose

setGeneric("transpose", function(tensor, foreignBackendOpt=FALSE) {
    result <- t(tensor)
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("transpose", c("MatrixLibAdapter", "ANY"), function(tensor, foreignBackendOpt=FALSE) {
    result <- transpose(tensor@matrix, FALSE)#tensor@foreignBackend)
    return(result);
})

# diagonal

setGeneric("diagonal", function(tensor, foreignBackendOpt=FALSE) {
    result <- Diagonal(x=as.vector(tensor))
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("diagonal", c("MatrixLibAdapter", "ANY"), function(tensor, foreignBackendOpt=FALSE) {
    result <- diagonal(tensor@matrix, FALSE); #tensor@foreignBackend)
    return(result);
})

# crossProduct

setGeneric("crossProduct", function(tensor, foreignBackendOpt=FALSE) {
    result <- crossprod(tensor)
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("crossProduct", c("MatrixLibAdapter", "ANY"), function(tensor, foreignBackendOpt=FALSE) {
    result <- crossProduct(tensor@matrix, tensor@foreignBackend)
    return(result);
})

# elementWiseSqrt

setGeneric("elementWiseSqrt", function(tensor, foreignBackendOpt=FALSE) {
    result <- (tensor)^{1/2}
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("elementWiseSqrt", c("MatrixLibAdapter", "ANY"), function(tensor, foreignBackendOpt=FALSE) {
    result <- elementWiseSqrt(tensor@matrix, FALSE);#tensor@foreignBackend)
    return(result);
})

# crossProductDuo 

setGeneric("crossProductDuo", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- crossprod(tensor, otherMatrix)
    #newTensor <- MatrixLibAdapter(matrix=result, foreignBackend=foreignBackendOpt)
    #return(newTensor);
    return(result);
})

setMethod("crossProductDuo", c("MatrixLibAdapter", "MatrixLibAdapter", "ANY"), function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- crossProductDuo(tensor@matrix, otherMatrix@matrix, FALSE);#tensor@foreignBackend)
    return(result);
})

setMethod("crossProductDuo", c("MatrixLibAdapter", "ANY", "ANY"), function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- crossProductDuo(tensor@matrix, otherMatrix, FALSE);#tensor@foreignBackend)
    return(result);
})


# Splice

setGeneric("splice", function(tensor, rowBeg, rowEnd, colBeg, colEnd) {

    rowBeg <- rowBeg + 1
    rowEnd <- rowEnd + 1
    colBeg <- colBeg + 1
    colEnd <- colEnd + 1

    result <- tensor[rowBeg:rowEnd, colBeg:colEnd]
    #newTensor <- MatrixLibAdapter(result)
    #return(newTensor)
    return(result);
})
 
setMethod("splice", 
    c("MatrixLibAdapter", "numeric", "numeric", "numeric", "numeric"),
    function(tensor, rowBeg, rowEnd, colBeg, colEnd) {
        result <- splice(tensor@matrix, rowBeg, rowEnd, colBeg, colEnd)
        #result@foreignBackend <- tensor@foreignBackend
        return(result);
})

setGeneric("removeAbstractions", function(tensor){
    if(tensor@foreignBackend){
        return(as.matrix(tensor@matrix))
    }
    
    val = tensor@matrix
    return(val);
})

# =============== The following needs to be abstracted away
setMethod("^", c("MatrixLibAdapter", "ANY"), function(e1, e2) {
    result <- e1@matrix ^ e2
    return(result)
})

setMethod("^", c("ANY", "MatrixLibAdapter"), function(e1, e2) {
    result <- e1 ^ e2@matrix
    return(result)
})

setMethod("%*%", c("MatrixLibAdapter", "ANY"), function(x, y) {
    result <- x@matrix %*% y
    return(result)
})

setMethod("%*%", c("ANY", "MatrixLibAdapter"), function(x, y) {
    result <- x %*% y@matrix
    return(result)
})

setMethod("*", c("MatrixLibAdapter", "ANY"), function(e1, e2) {
    result <- e1@matrix *  e2
    return(result)
})

setMethod("*", c("ANY", "MatrixLibAdapter"), function(e1, e2) {
    result <- e1 * e2@matrix
    return(result)
})

# ============================= Minor ops
setGeneric("invert", function(e1){
    result <- MatrixLibAdapter(matrix=(1 / e1@matrix))
    return(result)
})

setGeneric("divisionArr", function(e1, e2){
    result <- e2 / e1
    result <- MatrixLibAdapter(matrix=result)
    return(result)
})

setMethod("divisionArr", c("MatrixLibAdapter", "ANY"), function(e1, e2){
    return(divisionArr(e1@matrix, e2))
})

setMethod("divisionArr", c("MatrixLibAdapter", "MatrixLibAdapter"), function(e1, e2){
    return(divisionArr(e1@matrix, e2@matrix))
})
