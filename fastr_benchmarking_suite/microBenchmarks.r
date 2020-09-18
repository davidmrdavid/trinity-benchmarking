doScalarAddition <- function(testMatrix) {
    result <- testMatrix + 42
    return(result)
}

doScalarMultiplication <- function(testMatrix) {
    result <- testMatrix * 42
    return(result)
}

doElementWiseSum <- function(testMatrix) {
    result <- sum(testMatrix)
    return(result)
}

doRowWiseSum <- function(testMatrix) {
    result <- rowSums(testMatrix)
    return(result)
}

doColumnWiseSum <- function(testMatrix) {
    result <- colSums(testMatrix)
    return(result)
}

doLeftMatrixMultiplication <- function(testMatrix, otherMatrix) {
    result <- testMatrix %*% otherMatrix
    return(result)
}

doTransLeftMatrixMultiplication <- function(testMatrix, otherMatrix) {
    result <- t(testMatrix) %*% otherMatrix
    return(result)
}

doRightMatrixMultiplication <- function(testMatrix, otherMatrix) {
  result <- otherMatrix %*% testMatrix
  return(result)
}

doCrossProduct <- function(testMatrix) {
  result <- crossprod(testMatrix);
  return(result);
}
