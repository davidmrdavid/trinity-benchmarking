library(matrixStats)
source("morpheusR.r")

doLogisticRegression <- function(Data, max_iter, winit, gamma0, Target) {
    w = winit;
    for( k in 1:max_iter)
    {
        w =  gamma0 * (t(Data) %*%  (Target/(1+(2.78 ^ (Data%*%w))))); 
    }
    return(list(w))
}


doLinearRegression <- function(Data, Max_Iter, winit, gamma0, Target) {
    w = winit;
    for( k in 1:Max_Iter )
    {
        w = w - gamma0 * (t(Data) %*% ( Data %*% w - Target ));
    }
    return(list(w));
}


doKMeansClustering <- function(Data, Max_Iter, Center_Number, k_center, nS) {

    All1 = matrix(1, nS, 1);
    All1_k = matrix(1,1,Center_Number);
    All1_C = t(matrix(1,1,nrow(k_center)));	
    T2 = rowSums(Data^{2}) %*% All1_k;
    T22 = Data * 2;
    for( k in 1: Max_Iter )
    {
        Dist = T2 - as.matrix(T22 %*% k_center ) +  All1 %*% colSums(k_center ^2);
        YA = (Matrix(Dist == (rowMins((Dist)) %*% All1_k),sparse=TRUE))+0;
        arg1 <- t(Data) %*% YA;
        arg2 <- All1_C %*% colSums(YA);
        k_center = as.matrix(  ( t(Data) %*% YA ) / ( (All1_C) %*% colSums(YA) )  );
    }
    return(list(k_center ,YA));
}



doGNMFClustering <- function(X, Max_Iter, winit, hinit){
  w <- winit;
  h <- hinit;
  for(k in 1:Max_Iter) { 
    numerator <- (t(w) %*% X);
    denominator <- crossprod(w) %*% h;
    arg <- numerator / denominator;
    h <- (h * arg);
    # Updating w
    w <- (w * ((X %*% t(h)) / (w %*% (h %*% t(h)))));
  }
  return(list(h, w))
}


