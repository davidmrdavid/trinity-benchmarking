library(Matrix)

MatrixLibAdapter2 <- setClass(
  "MatrixLibAdapter2",
   slot = c(
     Sparse = "logical",
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
     scalarMultiplication = "ANY",
     scalarExponentiation = "ANY",
     elementWiseSum = "ANY"
   )
)


setMethod("initialize", "MatrixLibAdapter2",

    function(.Object, Sparse) {
    .Object@Sparse = Sparse
	.Object@rowSum = function(x) {
            return(rowSum(x));
        }
        .Object@columnSum = function(x) {
            return(columnSum(x));
        }
        .Object@leftMatrixMultiplication = function(x, y) {
            tStart <- as.numeric(Sys.time())*1000;
            z <- leftMatrixMultiplication(x, y)
            
            if(!Sparse & (toString(class(x)) == "ngCMatrix" | toString(class(y)) == "ngCMatrix")){ #SPARSE check
                z <- as.matrix(z)
            }
            dims <- dim(z)
            if(dims[2] == 1){
                z <- as.numeric(z)
            }
            tEnd <- as.numeric(Sys.time()) * 1000;
            time <- tEnd - tStart;
            print(sprintf("LMM: %i", time))
            return(z)
        }
        .Object@scalarAddition = function(x, y){
            return(scalarAddition(x,y));
  
        }
        .Object@scalarMultiplication = function(x, y){
            return(scalarMultiplication(x,y));
  
        }
        .Object@scalarExponentiation = function(x, y){
            return(scalarExponentiation(x, y))
        }
        .Object@rightMatrixMultiplication = function(x, y) {
            print(object.size(x), units = "auto");
            print(object.size(y), units = "auto");
            tStart <- as.numeric(Sys.time())*1000;
            z <- rightMatrixMultiplication(x, y)
            
            if(!Sparse & (toString(class(x)) == "ngCMatrix" | toString(class(y)) == "ngCMatrix")){ #SPARSE check
                z <- as.matrix(z)
            }
            dims <- dim(z)
            if(dims[2] == 1){
                z <- as.numeric(z)
            }
            tEnd <- as.numeric(Sys.time()) * 1000;
            time <- tEnd - tStart;
            print(sprintf("RMM: %i", time))
            return(z)
        }
        .Object@columnWiseAppend = function(x, y) {
            return(columnWiseAppend(x, y));
        }
        .Object@transpose = function(x) {
            return(transpose(x));
        }
        .Object@crossProduct = function(x) {
            z <- crossProduct(x)
            if(!Sparse  & toString(class(x)) == "ngCMatrix"){
                z <- as.matrix(z)
            }
            return(z);
        }
        .Object@crossProductDuo = function(x, y) {
            z <- crossProductDuo(x, y)
            if(!Sparse & (toString(class(x)) == "ngCMatrix" | toString(class(y)) == "ngCMatrix")){
                z <- as.matrix(z)
            }
            return(z);
        }
        .Object@elementWiseSum = function(x) {
            return(elementWiseSum(x));
        }
	.Object@elementWiseSqrt = function(x) {
            return(elementWiseSqrt(x));
        }
        .Object@rowWiseAppend = function(x, y) {
            return(rowWiseAppend(x, y));
        }
        .Object@diagonal = function(x) {
            return(diagonal(x));
        }
        .Object@splice = function(x, rowBeg, rowEnd, colBeg, colEnd) {
            return(splice(x, rowBeg, rowEnd, colBeg, colEnd));
        }

        .Object@matrixAddition = function(x, y) {
            cast = FALSE;
            
            if((class(x) != "numeric")){
                d <- dim(x);
                x = as.numeric(x);
                cast = TRUE;
            }
            if((class(y) != "numeric")){
                d <- dim(y);
                y = as.numeric(y);
                cast = TRUE;
            }
            
            z <- matrixAddition(x, y);
            if(cast){
                dim(z) <- d;
            }
            return(z);
        }
        .Object@getNumCols = function(x) {
            return(getNumCols(x));

        }
        .Object@getNumRows = function(x) {
            z <- getNumRows(x)
            return(z);

        }
        return(.Object)
})



# MatrixLibAdapter methods ====================================================


setGeneric("rightMatrixMultiplication", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {    
    return(tensor %*% otherMatrix);
})


setGeneric("leftMatrixMultiplication", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    return(otherMatrix %*% tensor);
})

setGeneric("scalarAddition", function(tensor, otherMatrix) {
    return(tensor + otherMatrix);
})

setGeneric("scalarExponentiation", function(tensor, number) {
        return(tensor ^ number); #number ^ tensor, would fail (?)
})

setGeneric("scalarMultiplication", function(tensor, scalar) {
    return(tensor * scalar);
})


setGeneric("crossProduct", function(tensor) {
    return(tensor * normMatrix);
})


setGeneric("rowSum", function(tensor, foreignBackendOpt=FALSE) {
    return(rowSums(tensor));
})

setGeneric("columnSum", function(tensor, foreignBackendOpt=FALSE) {
    return(colSums(tensor));
})

setGeneric("elementWiseSum", function(tensor) {
    return(sum(tensor));
})


setGeneric("rowWiseAppend", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    return(rbind(tensor, otherMatrix));
})

setGeneric("columnWiseAppend", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    return(cbind(tensor, otherMatrix));
})

setGeneric("matrixAddition", function(tensor, otherMatrix) {
    return(tensor + otherMatrix);
})

setGeneric("getNumRows", function(tensor) {
    return(nrow(tensor));
})

setMethod("getNumRows", c("MatrixLibAdapter"), function(tensor) {
    result <- getNumRows(tensor@matrix)
    return(result);
})

setGeneric("getNumCols", function(tensor) {
    return(ncol(tensor));
})

setMethod("getNumCols", c("MatrixLibAdapter"), function(tensor) {
    return(getNumCols(tensor@matrix));
})

setGeneric("transpose", function(tensor, foreignBackendOpt=FALSE) {
    return(t(tensor));
})

setGeneric("diagonal", function(tensor, foreignBackendOpt=FALSE) {
    # TODO: why do we need the `as.vector`
    return(Diagonal(x=as.vector(tensor)));
})


setGeneric("crossProduct", function(tensor, foreignBackendOpt=FALSE) {
    return(crossprod(tensor));
})


setGeneric("elementWiseSqrt", function(tensor, foreignBackendOpt=FALSE) {
    return((tensor)^{1/2});
})


setGeneric("crossProductDuo", function(tensor, otherMatrix, foreignBackendOpt=FALSE) {
    result <- crossprod(tensor, otherMatrix)
    return(result);
})


setGeneric("splice", function(tensor, rowBeg, rowEnd, colBeg, colEnd) {

    rowBeg <- rowBeg + 1
    rowEnd <- rowEnd + 1
    colBeg <- colBeg + 1
    colEnd <- colEnd + 1

    return(tensor[rowBeg:rowEnd, colBeg:colEnd]);
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
