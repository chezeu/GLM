 
########### logistic equations 
# context 1

# EM algorithm for logistic regression with probabilistic linkage (Context 1)
 
# logistic probability 
P_matrix <- function( X12,X34, p1,p2, nA, nB, beta) {
  
#  p = p1+p2
 # X1=as.matrix(X1,ncol=p1,nrow=nB)
  #X12=as.matrix(X12,ncol=p2,nrow=nA)
  
  eta1 <- X34 %*% beta[-c(1:p1)]        # length nB
  eta2 <- X12 %*% beta[1:p1]    # length nA
  # form matrix eta_mat of size nA x nB: eta_mat[i,j] = eta1[i] + eta2[j]
  eta_mat <- t(matrix(rep(as.numeric(eta1), times = nA), nrow = nB)) +
    t( matrix(rep(as.numeric(eta2), each  = nB), nrow = nB))
  P_mat <- 1 / (1 + exp(-eta_mat))   # P_{ij}
  
  return(P_mat)
}

#posterior probability
post_proba <- function( Z, P_mat, Q){

  #nA=nrow(Q)
 # nB=ncol(Q)
 
  Z_mat <- matrix(Z, nrow = nA, ncol = nB, byrow = TRUE)
  numerator <- (P_mat^Z_mat) * ((1 - P_mat)^(1 - Z_mat)) * Q # nA x nB
 # denominator <- apply(numerator,1, sum)
  denominator <- rowSums(numerator)
  if (denominator == 0 || !is.finite(denominator)) denominator <- .Machine$double.eps
  
  matrix_proba_Z <- numerator/denominator # numerators in log scale
  
  return(matrix_proba_Z)
}

####### estimating equation
logistic_equation_1 <- function( beta,p1,p2, matrix_proba_Z, Z, X12,X34){

  #p = p1+p2
  #X3=as.matrix(X3,ncol=p1,nrow=nB)
  #X12=as.matrix(X12,ncol=p2,nrow=nA)
  #Z= as.matrix(Z, ncol=1,nrow= nB)
  # nA= nrow(matrix_proba)
  #nB= ncol(matrix_proba)
  
  P_mat <- P_matrix( X12,X34,p1,p2, nA, nB, beta)
  
  Z_mat <- matrix(Z, nrow = nrow(P_mat), ncol = length(Z), byrow = TRUE)
  W <- matrix_proba_Z * (Z_mat - P_mat)
  
  W <- matrix_proba_Z * (Z_mat - P_mat)   # nA x nB
  score2 <- colSums(W %*% X34) 
  score1 <- colSums(t(W) %*% X12)   # ??2
  
 s=  c(score1, score2)
  return(s)
}

# ##########solving equation

logistic_1 <- function(beta_ini, Z, X12,X34, matrix_proba_Z, maxiter = 20){
  
  
  f <- function(x){
    logistic_equation_1( beta=x ,p1,p2, matrix_proba_Z, Z, X12,X34)
  }
  fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list(coef = beta0, converge = converge, iterations = iterations))
}


### iterrations ###########################################

logistic_itteration<-function(beta0, p1,p2, Z, X12,X34, Q, tol= 1e-6, maxits = 300){
  
  X34=as.matrix(X34)
  X12=as.matrix(X12)
  nA=nrow(X12)
  nB=nrow(X34)
  p = p1+p2
  Z= as.matrix(Z )

  beta_ini = rep(0, p)
  it = 0
  converge = FALSE
  
  if(nrow(Z) != nB) stop("length(Z) must equal nrow(X12)")
  if(!all(dim(Q) == c(nA, nB))) stop("Q must be nA x nB")
  
  
  while ((!converge) && (it < maxits)){ 
    
    beta0_old = beta0
    
    #logistic proba
    P_mat = P_matrix( X12,X34,p1,p2, nA, nB, beta0_old)
    
    #expectation  
    matrix_proba_Z = post_proba( Z, P_mat, Q)
    
    #maximization
    estime = logistic_1(beta_ini=beta0_old, Z, X12, X34,matrix_proba_Z, maxiter = 20)
    beta0 = estime$coef
    
    converge <- sqrt(sum((beta0 - beta0_old)^2)) < tol
    
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
  return(list(beta0=beta0, converge= converge, it=it))
}
 ########### naive estimation
      
##true estimation
 
  