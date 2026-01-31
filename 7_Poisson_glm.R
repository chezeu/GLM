
# EM algorithm for logistic regression with probabilistic linkage (Context 1)
 
########
### proba aposteriories (pi)

Func_prob_V<-function(beta,V, X12,X34,p1,p2, Q){
  
  #X3=as.matrix(X3,ncol=p2,nrow=nB)
  # X12=as.matrix(X12,ncol=p1,nrow=nA)
  # nA=nrow(X12)
  # nB=nrow(X3)
  # p = p1+p2
  
  eta1 <- X34 %*% beta[-c(1:p1)]        # length nB
  eta2 <- X12 %*% beta[1:p1]    # length nA
  # form matrix eta_mat of size nA x nB: eta_mat[i,j] = eta1[i] + eta2[j]
  eta_mat <- t(matrix(rep(as.numeric(eta1), times = nA), nrow = nB)) +
    t( matrix(rep(as.numeric(eta2), each  = nB), nrow = nB))
  
  v_matrix <- matrix(V, nrow = nA, ncol = length(V), byrow = TRUE)
 facto_V = matrix( lgamma(V + 1), nrow = nA, ncol = length(V), byrow = TRUE) # log(V!) computed using lgamma for numerical stability
 log_pois <- v_matrix * eta_mat - exp(eta_mat) - facto_V
  
 logQ <- matrix(-Inf, nA, nB)
 logQ[Q > 0] <- log(Q[Q > 0]) # avoid log(0)
 log_num <- logQ + log_pois # in order to do exp(log(num))
 
 #matrix_proba_V <- matrix(0, nA, nB)
 
# for (i in 1:nA) {
 #  m <- max(log_num[i, ], finite = TRUE)
   
  # if (!is.finite(m)) next  # ligne impossible
   
   #log_denom <- m + log(sum(exp(log_num[i, ] - m)))
  # matrix_proba_V[i, ] <- exp(log_num[i, ] - log_denom)
 #}
 

 
 # stabilisation softmax each ligne  
 
 log_num <- log_num - apply(log_num, 1, max)
 
 num <- exp(log_num)
 denom <- rowSums(num)
 denom[denom == 0] <- .Machine$double.eps
 matrix_proba_V = num /denom
 
 return(matrix_proba_V)
 #num= Q * exp(log_pois)
 #denominateur = rowSums (num )
#  if (denominateur == 0 || !is.finite(denominateur)) denominateur <- .Machine$double.eps
#  matrix_proba_V = num /denominateur
  
} 
 ###################
#  estimating equation ########################
Poisson_equation_1 <- function(beta,matrix_proba_V,V, 
                           X12,X34,p1,p2) {
   
  #X34=as.matrix(X34,ncol=p2,nrow=nB)
  # X12=as.matrix(X12,ncol=p1,nrow=nA)
  # nA=nrow(X12)
  # nB=nrow(X3)
  # p = p1+p2
  
  eta1 <- X34 %*% beta[-c(1:p1)]        # length nB
  eta2 <- X12 %*% beta[1:p1]    # length nA
  # form matrix eta_mat of size nA x nB: eta_mat[i,j] = eta1[i] + eta2[j]
  eta_mat <- t(matrix(rep(as.numeric(eta1), times = nA), nrow = nB)) +
    t( matrix(rep(as.numeric(eta2), each  = nB), nrow = nB))
  eta_mat <- pmin(eta_mat, 20)# enough for the poisson process
  
  ####### 
  v_matrix <- matrix(V, nrow = nA, ncol = length(V), byrow = TRUE)

  vi_beta= v_matrix- exp(eta_mat)  # matrix nA*nB
   
  H <- rep(0, p)
  
  for (i in 1:nA) {
    for (j in 1:nB) {
      
      X_ij = c (X12[i,],X34[j,]  )
      
      H <- H +  matrix_proba_V [i, j] * X_ij * vi_beta[i,j]
    }
  }
  
  return(as.numeric(H))
  
}

#####################
# ###########solve the equation 

Poisson_1<- function( matrix_proba_V, V, X12,X34,
                  p1,p2,beta_ini,maxiter = 20){
  
  f <- function(x){
    val<-Poisson_equation_1( beta=x,matrix_proba_V,V, X12,X34,p1,p2)
    
    if (!is.numeric(val) || length(val) != length(x)) {
      stop(" Poisson_equation_1 must return a numeric vector of same length as x")
    }
    return(val)
  }
  fit_manual <- nleqslv( beta_ini,f, method = "Newton")
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1) )
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}
################

### iterrations ###########################################

#valeurs initials

Poisson_itteration<-function(beta0,V, X12,X34,
                         p1,p2, Q,tol= 1e-6, maxits = 300){
  
  X34=as.matrix(X34)
  X12=as.matrix(X12)
  V= as.matrix(V)
  nA=nrow(X12)
  nB=nrow(X34)
  p = p1+p2
  
  beta_ini=rep(0,p)
  it = 0
  converge = FALSE
  
  while ((!converge) && (it < maxits)){ 
     
    beta0_old = beta0
    
    #expectation  
      
    matrix_proba_V= Func_prob_V (beta0_old ,V, X12,X34,p1,p2, Q)
    
    #maximization
    estime = Poisson_1 ( matrix_proba_V,V, X12,X34,p1,p2,beta_ini,maxiter = 20)
    beta0 = estime$beta0
    beta_ini = beta0_old 
    
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
  return(list(beta0=beta0 , matrix_proba_V=matrix_proba_V, 
              converge= converge, it=it))
}
 

