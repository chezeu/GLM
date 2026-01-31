
####expectation step
#posterior probability


func_proba_Y <- function( Y, X12,X34,p1,p2, variance1, beta, Q){
  
  
  eta_mat= expo_beta  (X12,X34, p1,p2 )  
  Y_mat <- matrix(Y, nrow = nA, ncol = nB, byrow = TRUE)
  #
  #logQ <- matrix(-Inf, nA, nB)
  #logQ[Q > 0] <- log(Q[Q > 0]) # avoid log(0)
  #  log_num <- -(1/2*variance1)*(Y_mat - eta_mat)^2 + logQ
  
  # stabilisation softmax each ligne  
  #log_num <- log_num - apply(log_num, 1, max)
  #num <- exp(log_num)
  
  num <- Q * exp( -(1/2*variance1)*(Y_mat - eta_mat)^2)
  
  denom <- rowSums(num)
  denom[denom == 0] <- .Machine$double.eps
  
  matrix_proba_Y=  num / denom
  
  return(matrix_proba_Y)
}

###########
linear_equation_1 <- function(p1,p2, Y,Q,matrix_proba_Y,
                              X12,X34){
  nA <- nrow(X12)
  nB <- nrow(X34)
  p  <- length(p1)  +length(p2) 
  
  # normal equations A*beta=b
  A <- matrix(0, p, p)
  b <- numeric(p)
  
  for (i in 1:nA) {
    W <- diag(matrix_proba_Y[i, ])
    
    Xi <- cbind(matrix(X12[i, ], nB, p1, byrow = TRUE),X34  )
    
    A <- A + t(Xi) %*% W %*% Xi
    b <- b + t(Xi) %*% (matrix_proba_Y[i, ] * Y)
  }
  
  if (rcond(A) < 1e-12) stop("Normal equation matrix is singular")
  
  beta_new <- solve(A, b)
  
  # sigma² 
  eta <- X34 %*% beta_new[p2]  # size nB
  sigma2 <- 0
  
  for (i in 1:nA) {
    mu <- eta + drop(X12[i, ] %*% beta_new[p1])
    sigma2 <- sigma2 + sum(matrix_proba_Y[i, ] * (Y - mu)^2)
  }
  
  sigma2 <- sigma2 / nA
  
  return(list(beta = beta_new, sigma2=sigma2))
}

### iterrations ###########################################

linear_itteration<-function(beta0,variance1, p1,p2, Y, X12,X34, Q, 
                            tol= 1e-6, maxits = 300){
  
  X34=as.matrix(X34)
  X12=as.matrix(X12)
  Y= as.matrix(Y)
  nA=nrow(X12)
  nB=nrow(X34)
  p = length(p1)+length(p2)
  it = 0
  converge = FALSE
  
  if(nrow(Y) != nB) stop("length(Y) must equal nrow(X34)")
  if(!all(dim(Q) == c(nA, nB))) stop("Q must be nA x nB")
  
  
  while ((!converge) && (it < maxits)){ 
    
    beta0_old = beta0
    variance1_old=variance1
    
    #expectation  
    matrix_proba_Y = func_proba_Y( Y, X12,X34,p1,p2, 
                                   variance1_old, beta=beta0_old, Q)
    
    #maximization
    estime = linear_equation_1 ( p1,p2, Y,Q,matrix_proba_Y,
                                 X12,X34)
    beta0 = estime$beta
    variance1 =estime$sigma2
    
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
  return(list(beta0=beta0, variance1= variance1,
              converge= converge, it=it))
} 

