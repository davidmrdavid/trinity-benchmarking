source("matrixLibAdapter.r")

# Adapter, stores NM as member variable
NormalizedMatrix <- setClass(
   "NormalizedMatrix",
   slot = c(morpheus = "ANY")
)

# Adapter constructor
asNormalizedMatrix <- function(S, Ks, Rs, Sparse=FALSE) {

    # Obtain NM constructor, execute it, store it in adapter object,
    # return adapter
    morpheusBuilder <- eval.polyglot("morpheusDSL", "")
    # TODO: check is S is empty
    Sempty <- FALSE
    if(nrow(S)*ncol(S) == 0){
        Sempty <- TRUE
    }
    avatar <- MatrixLibAdapter2(Sparse=Sparse)
    morpheus <- morpheusBuilder@build(S, Ks, Rs, Sempty, avatar)
    normMatrix <- NormalizedMatrix(morpheus=morpheus)
    result <- normMatrix.scalarAddition(40)
    return(normMatrix)
}

setMethod("show", "NormalizedMatrix", function(object){
    print("A NormalizedMatrix")
})

setMethod("+", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    result <- e2@morpheus@scalarAddition(e1)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("+", c("NormalizedMatrix", "numeric"), function(e1, e2) {
    result <- e1@morpheus@scalarAddition(e2)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("-", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    result <- e2@morpheus@scalarAddition(e1)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("-", c("NormalizedMatrix", "numeric"), function(e1, e2) {
    result <- e1@morpheus@scalarAddition(e2)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("*", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    result <- e2@morpheus@scalarMultiplication(e1)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("*", c("NormalizedMatrix", "numeric"), function(e1, e2) {
    result <- e1@morpheus@scalarMultiplication(e2)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("/", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    preppedArg <- 1/e2
    e1@morpheus = e1@morpheus@scalarMultiplication(preppedArg)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("/", c("NormalizedMatrix", "numeric"), function(e1, e2) {
    preppedArg <- 1/e1
    result <- e2@morpheus@scalarMultiplication(preppedArg)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

# TODO: this is wrong!
setMethod("^", c("NormalizedMatrix", "numeric"), function(e1, e2) {
    result <- e1@morpheus@scalarExponentiation(e2)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("^", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    result <- e2@morpheus@scalarExponentiation(e1)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("%*%", c("ANY", "NormalizedMatrix"), function(x, y) {
    result <- y@morpheus@rightMatrixMultiplication(x)
    return(result);
})


setMethod("%*%", c("NormalizedMatrix", "ANY"), function(x, y) {
    result <- x@morpheus@leftMatrixMultiplication(y)
    return(result);
})


setMethod("sum", c("NormalizedMatrix"), function(x) {
    result <- x@morpheus@elementWiseSum()
    return(result);
})

setMethod("rowSums", c("NormalizedMatrix"), function(x) {
    result <- x@morpheus@rowSum()
    return(result)
})

setMethod("colSums", c("NormalizedMatrix"), function(x) {
    result <- x@morpheus@columnSum()
    return(result)
})


't.NormalizedMatrix' <- function(x) {
    result <- x@morpheus@transpose()
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
}


setMethod("crossprod", c("NormalizedMatrix", "ANY"), function(x, y = NULL) {
    result <- x@morpheus@crossProduct();
    return(result);
})
