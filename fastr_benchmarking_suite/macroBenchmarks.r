library(matrixStats)
source("morpheusR.r")

doLogisticRegression <- function(Data, max_iter, winit, gamma0, Target) {
    w = winit;
    Max_Iter = 2;
    for( k in 1:max_iter)
    {
        tStart <- as.numeric(Sys.time())*1000;
        w =  gamma0 * (t(Data) %*%  (Target/(1+(2.78 ^ (Data%*%w))))); 
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("1: %d", time));
    }
    return(list(w))
}


doLinearRegression <- function(Data, Max_Iter, winit, gamma0, Target) {
    w = winit;
    Max_Iter = 1;
    for( k in 1:Max_Iter )
    {      
        tStart <- as.numeric(Sys.time())*1000;
        d = t(Data);
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("1: %d", time));
        
        #print(toString(class(Data))); #Normal
        #print(toString(class(w))); #dge
        
        tStart <- as.numeric(Sys.time())*1000;
        d2 =Data %*% w; 
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("2: %d", time));
        
        #print(toString(class(d2))); #Numeric Dge
        #print(toString(class(Target))); #Matrix
        
        tStart <- as.numeric(Sys.time())*1000;
        d2 = d2 - Target;
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("3: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        d = d %*% d2;
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("4: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        d = gamma0 * d
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("5: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        w = w - d
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("6: %d", time));
        
    }
    return(list(w));
}


doKMeansClustering <- function(Data, Max_Iter, Center_Number, k_center, nS) {

    #print(class(Data))
    #X
    #print("20")
    All1 = matrix(1, nS, 1);
    All1_k = matrix(1,1,Center_Number);
    All1_C = t(matrix(1,1,nrow(k_center)));	
    T2 = rowSums(Data^{2}) %*% All1_k;
    print("T2");
    print(class(T2));
    #print(W)
    T22 = Data * 2;
    Max_Iter = 1;
    for( k in 1: Max_Iter )
    {
        #print(class(Data))
        #print(class(Data^{2}))
        #print(class(rowSums(Data^{2})))
        #print(class(All1_k))
        #print(class(T2))
        
        #b =  colSums(k_center ^2)#(T22 %*% k_center )
        #print("###############")
        #print(class(b))
        #print(X)
        tStart <- as.numeric(Sys.time())*1000;
#         print(dim());
#         print(dim(k_center));
        
        tmp = T22 %*% k_center;
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("1: %d", time));
        
        
        tStart <- as.numeric(Sys.time())*1000;
        Dist = T2 - as.matrix(tmp);
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("2A: %d", time));
        #print(class(Dist))
        #print(X)
        tStart <- as.numeric(Sys.time())*1000;
        Dist = Dist +  All1 %*% colSums(k_center ^2);
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("2B: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        YA = (Matrix(Dist == (rowMins((Dist)) %*% All1_k),sparse=TRUE))+0;
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("3: %d", time));
        
        #print(Z)
        tStart <- as.numeric(Sys.time())*1000;
        arg1 <- t(Data) %*% YA;
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("4: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        arg2 <- All1_C %*% colSums(YA);
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("5: %d", time));
        
        tStart <- as.numeric(Sys.time())*1000;
        k_center = as.matrix(  ( t(Data) %*% YA ) / ( (All1_C) %*% colSums(YA) )  );
        tEnd <- as.numeric(Sys.time()) * 1000;
        time <- tEnd - tStart;
        print(sprintf("6: %d", time));
        
        
    }
    return(list(k_center ,YA));
}



doGNMFClustering <- function(X, Max_Iter, winit, hinit){
  
  w <- winit;
  h <- hinit;
  print(dim(w));
  print(dim(h));
  Max_Iter = 10;
  for(k in 1:Max_Iter) { 
#     numerator <- (t(w) %*% X);
#     denominator <- crossprod(w) %*% h;
#     arg <- numerator / denominator;
#     h <- (h * arg);
#     # Updating w
#     w <- (w * ((X %*% t(h)) / (w %*% (h %*% t(h)))));
    tStart <- as.numeric(Sys.time())*1000;
    numerator <- (t(w) %*% X);
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("1: %d", time));
    
    tStart <- as.numeric(Sys.time())*1000;
    denominator <- crossprod(w) %*% h;
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("2: %d", time));
    
    tStart <- as.numeric(Sys.time())*1000;
    arg <- numerator / denominator;
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("3: %d", time));
    
    tStart <- as.numeric(Sys.time())*1000;
    h <- (h * arg);
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("4: %d", time));
    
#     print(toString(class(X)));
#     print(toString(class(h)));
    tStart <- as.numeric(Sys.time())*1000;  
    tmp1 = X %*% t(as.matrix(h));
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("5: %d", time));
     
    tStart <- as.numeric(Sys.time())*1000;
    tmp2 = h %*% t(h);
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("6: %d", time));
    
    tStart <- as.numeric(Sys.time())*1000;
    tmp2 = w %*% tmp2;
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("7: %d", time));
    
    tStart <- as.numeric(Sys.time())*1000;
    w <- w * tmp1;
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("8: %d", time));
            
    tStart <- as.numeric(Sys.time())*1000;
    w <- w / tmp2;  
    tEnd <- as.numeric(Sys.time()) * 1000;
    time <- tEnd - tStart;
    print(sprintf("9: %d", time));
    # Updating w
    #w <- (w * ((X %*% t(h)) / (w %*% (h %*% t(h)))));
  }
  return(list(h, w))
}


