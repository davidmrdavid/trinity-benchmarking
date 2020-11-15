source("matrixLibAdapter.r")

# Adapter, stores NM as member variable
NormalizedMatrix <- setClass(
   "NormalizedMatrix",
   slot = c(morpheus = "ANY")
)

# Adapter constructor
asNormalizedMatrix <- function(S, Ks, Rs) {

    # Wrap tensors in adapters 
    tensorS <- S #MatrixLibAdapter(matrix=S)
    tensorKs <- Ks #lapply(Ks, function(x){MatrixLibAdapter(matrix=x)}) 
    tensorRs <- Rs #lapply(Rs, function(x){MatrixLibAdapter(matrix=x)})

    # Obtain NM constructor, execute it, store it in adapter object,
    # return adapter
    morpheusBuilder <- eval.polyglot("morpheusDSL", "")
    # TODO: check is S is empty
    Sempty <- FALSE
    if(nrow(S)*ncol(S) == 0){
        Sempty <- TRUE
    }
    avatar <- MatrixLibAdapter2()
    morpheus <- morpheusBuilder@build(tensorS, tensorKs, tensorRs, Sempty, avatar)
    normMatrix <- NormalizedMatrix(morpheus=morpheus)
    return(normMatrix)
}

# I want the NormalizedMatrix to be hidden even further!!!!!!!
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
    print("AQUI")
    print(class(e2))
    print(class(e1))
    result <- e1@morpheus@scalarExponentiation(e2)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("^", c("numeric", "NormalizedMatrix"), function(e1, e2) {
    print("ALLA")
    result <- e2@morpheus@scalarExponentiation(e1)
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
})

setMethod("%*%", c("matrix", "NormalizedMatrix"), function(x, y) {
    tensorArg <- MatrixLibAdapter(x)
    result <- y@morpheus@rightMatrixMultiplication(tensorArg)
    return(removeAbstractions(result))
})

# TODO: add dgCMatrix everywhere ?

setMethod("%*%", c("NormalizedMatrix", "matrix"), function(x, y) {
    #tensorArg <- MatrixLibAdapter(y)
    result <- x@morpheus@leftMatrixMultiplication(y)
    #return(removeAbstractions(result))
    return(result);
})

setMethod("%*%", c("dgCMatrix", "NormalizedMatrix"), function(x, y) {
    #tensorArg <- MatrixLibAdapter(x)
    result <- y@morpheus@rightMatrixMultiplication(x)
    #return(removeAbstractions(result))
    return(result)
})

setMethod("%*%", c("NormalizedMatrix", "dgCMatrix"), function(x, y) {
    #tensorArg <- MatrixLibAdapter(y)
    result <- x@morpheus@leftMatrixMultiplication(y)
    #return(removeAbstractions(result))
    return(result);
})

setMethod("%*%", c("dgeMatrix", "NormalizedMatrix"), function(x, y) {
    #tensorArg <- MatrixLibAdapter(x)
    result <- y@morpheus@rightMatrixMultiplication(x)
    #return(removeAbstractions(result))
    return(result)
})

setMethod("%*%", c("NormalizedMatrix", "dgeMatrix"), function(x, y) {
    #tensorArg <- MatrixLibAdapter(y)
    result <- x@morpheus@leftMatrixMultiplication(y)
    #return(removeAbstractions(result))
    return(result);
})

setMethod("sum", c("NormalizedMatrix"), function(x) {
    result <- x@morpheus@elementWiseSum()
    return(as.vector(removeAbstractions(result)))
})

setMethod("rowSums", c("NormalizedMatrix"), function(x) {
    # Cast to vector, to match R's behaviour
    result <- x@morpheus@rowSum()
    return(as.vector(removeAbstractions(result)))
})

setMethod("colSums", c("NormalizedMatrix"), function(x) {
    result <- x@morpheus@columnSum()
    #return(as.vector(removeAbstractions(result)))
    return(result)
})


't.NormalizedMatrix' <- function(x) {
    result <- x@morpheus@transpose()
    newNormalizedMatrix <- NormalizedMatrix(morpheus=result)
    return(newNormalizedMatrix)
}


setMethod("crossprod", c("NormalizedMatrix", "ANY"), function(x, y = NULL) {
    result <- x@morpheus@crossProduct();
    #return(removeAbstractions(result));
    return(result);
})
