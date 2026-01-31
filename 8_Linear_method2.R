

#posterior probability

func_proba_Y_2 <- function( Y, X12,X34,p1,p2,
                            variance1, variance0, beta,mu2,matrix_id){
  
  #Y=as.matrix(Y)
  # X34=as.matrix(X34)
  #X12=as.matrix(X12)
  # nA =nrow(X12)
  #nB= nrow(X34)
  
  mu2X <- X34%*% mu2  # nB
  
  eta_mat= expo_beta  (X12,X34, p1,p2 )  
  
  Y_mat <- matrix(Y, nrow = nA, ncol = nB, byrow = TRUE)
  
  m = exp(-(1/2*variance0)*( Y-mu2X )^2) # nB
  m_mat = (1-matrix_id )* matrix(m, nrow = nA, ncol = nB, byrow = TRUE) # nA x nB
  
  num <- matrix_id * exp(-(1/2*variance1)*(Y_mat - eta_mat)^2)  #nA X nB
  #
  denom =  num  +  m_mat 
  
  denom[denom == 0] <- .Machine$double.eps
  
  matrix_proba_Y_2=  num / denom
  
  return(matrix_proba_Y_2)
}


linear_method_2 <- function( p1,p2, Y, matrix_proba_Y_2,
                             X12, X34){
  
  #nA =nrow(X12)
  # nB= nrow(X1)
  
  X_M <- matrix(0, nA * nB, ncol(X12) + ncol(X34)) #(nA*nB,p1+p2)
  idx <- 1
  
  for (i in 1:nA) {
    for (j in 1:nB) {
      X_M[idx, ] <- c(X12[i, ], X34[j, ] )
      idx <- idx + 1
    }
  }
  
  Y_M =NULL
  for (i in 1:nA) { 
    Y_M= rbind(Y_M,Y)
    i=i+1
  } 
  # vecteur diagonale des poids
  
  vect_pi = as.vector(t(matrix_proba_Y_2))
  
  ## v_j and V
  vj <- colSums(matrix_proba_Y_2) #nB
  V <- diag(vj)
  
  ## mu2
  mu2 <-as.vector( solve(
    t(X34) %*% V %*% X34,
    t(X34) %*% V %*% Y)
  )
  ## sigma0
  sigma_0 <- as.numeric( sum(vj * (Y - X34 %*% mu2)^2) / nA)
  
  #   beta
  
  beta=as.vector( solve( t(X_M) %*% (vect_pi * X_M), 
                         t(X_M) %*% (vect_pi * Y_M)) )
  
  #  sigma1
  res <- Y_M - X_M %*% beta
  sigma_1 <- as.numeric( (t(res) %*% ( vect_pi*res) )/ nA)
  
  return(list(beta =beta, mu2= mu2, sigma_1=sigma_1, sigma_0=sigma_0))
}

### iterrations ###########################################

Linear_itt_method2<-function(beta0,mu2_0,variance0,variance1, p1,p2, Y, X12, X34,  
                             matrix_id, tol= 1e-6, maxits = 300){
  
  X34=as.matrix(X34)
  X12=as.matrix(X12)
  Y= as.matrix(Y)
  nA=nrow(X12)
  nB=nrow(X34)
  p = length(p1)+length (p2)
  
  it = 0
  converge = FALSE
  
  if(nrow(Y) != nB) stop("length(Y) must equal nrow(X3)")
  if(!all(dim(matrix_id) == c(nA, nB))) stop("Q must be nA x nB")
  
  
  while ((!converge) && (it < maxits)){ 
    
    beta0_old = beta0
    mu0_old = mu2_0
    variance0_old=variance0
    variance1_old=variance1
    
    
    #expectation  
    matrix_proba_Y_2 = func_proba_Y_2( Y, X12,X34,p1,p2, variance1_old,
                                       variance0_old, beta0_old, mu0_old,matrix_id)
    
    #maximization
    estime = linear_method_2 ( p1,p2, Y,matrix_proba_Y_2,
                               X12,X34)
    beta0 = estime$beta
    variance0 =estime$sigma_0
    mu2_0 = estime$mu2
    variance1 =estime$sigma_1
    
    d_beta   <- sqrt(sum((beta0 - beta0_old)^2))
    d_variance1 <- sqrt(sum((variance1 - variance1_old)^2))
    converge <- (d_beta + d_variance1) < tol
    
    if (it == maxits && !converge) {
      cat("WARNING! NOT CONVERGENT WITH EM!", "\n")
      converge = FALSE
    }
    if (any(is.na(beta0))) {
      cat("WARNING! beta0 NOT AVAILABLE!", "\n")
      converge <- FALSE
    }
    
    it = it + 1
    
  }
  return(list(beta =beta0, mu2 =mu2_0 , variance0= variance0,  variance1= variance1,
              converge= converge, it=it))
} 








