
###############  matrix Q

compare_binary <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid( XB.k,XA.k)
    
    gamma.k = as.numeric(temp[,1]==temp[,2])
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  vect_gamma= vector()
  for (i in 1: nrow(comp_mat)) {
    if(sum(comp_mat [i,1:K])==K){vect_gamma[i]=1 
    }else {vect_gamma[i]=0 }
    
  }
  comp_mat [,K+1] = vect_gamma
 colnames(comp_mat)=c(paste("C_R", 1:K, sep=""), "vect_gamma" )
  
  return(comp_mat)
}

###############
Matrix_Q = function(K,datA,datB){
  
   nA= nrow(datA)
  nB=nrow(datB)
  idBA= datA[,K+1]
  
  comp_mat = compare_binary(datA, datB, K)
  comp_mat= data.frame(comp_mat) 
  vect_gamma= comp_mat$vect_gamma # comparison pairs
  
  matrix_id <- t( matrix(vect_gamma,ncol=nA,nrow=nB))
 # # common variables
  number = apply(matrix_id,1, FUN= function(x){length(which(x==1))}) # number of individual
  
  Q = matrix_id/number
  
  naive_id = apply(Q,1, which.max)
  
  return(list(Q=Q, naive_id=naive_id, matrix_id=matrix_id))
  
}
