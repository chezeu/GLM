
########### logistic equations 
# context 1

################### cumulative function

Funct_lambda2<-function(lambda0,Ts){
  #lambda0<- rep(0.1,length(event))
  nB = length(Ts)
  #lambda2 = vector()
  lambda2=numeric(nB)# vecteur nule
  
  for (j in 1:nB) {
    lambda2[j] =  sum(lambda0[which(Ts<=Ts[j])] )  
  }
  return(lambda2=lambda2)
}
########
### proba aposteriories (pi)

Func_prob_T<-function(beta,lambda2, event,X12,X34, p1,p2, Q){

  #X34=as.matrix(X34 )
 # X12=as.matrix(X12 )
 # nA=nrow(X12)
 # nB=nrow(X34)
 
  
  eta_mat= expo_beta  (X12,X34, p1,p2,beta )
   
  lambda_mat <-matrix(lambda2, nrow = nA, ncol = length(lambda2), byrow = TRUE)
 
  delta_ij <- matrix(event, nrow = nA, ncol = length( event), byrow = TRUE)
  
  # matrice pi(ij)
  numerat = ( exp(eta_mat) )^delta_ij * exp( -lambda_mat * exp(eta_mat) )*Q #size nA*nB
 
  denominateur = rowSums (numerat)
  denominateur[!is.finite(denominateur) | denominateur == 0] <- .Machine$double.eps
  matrix_proba_T = numerat/denominateur
  
  return( matrix_proba_T = matrix_proba_T)
}
#############
################ estimations ###################################
 
Function_lambda0<-function(matrix_proba_T,beta,Ts,event, X12,X34, p1,p2){
  
  #X34=as.matrix(X34,ncol=p2,nrow=nB)
 # X12=as.matrix(X12,ncol=p1,nrow=nA)
  #nA=nrow(X12)
 # nB=nrow(X34)
 # p = p1+p2 
  eta_mat= expo_beta  (X12,X34, p1,p2,beta )
  
  delta_ij <- matrix(event, nrow = nA, ncol = nB, byrow = TRUE)
  W = matrix_proba_T*delta_ij
  
    num <- numeric(nB)
    lambda0_denom <- numeric(nB)
  for (j in 1:nB) {
       t_d= Ts[j]
        n = which(Ts==t_d)
         num[j] = sum(W[,n]) 
  pi_exp_beta = matrix_proba_T * exp(eta_mat)  # matrix nA*nB
  
  risq = GetRiskSet (Ts[j], Ts)
  nrisk = length(risq)
  
   lambda0_denom[j]  = sum( pi_exp_beta [,risq] )# sum on j
                  
  }
    lambda0_denom[!is.finite(lambda0_denom) | lambda0_denom == 0] <- .Machine$double.eps
    
  lambda0_1 =  num / lambda0_denom
 # lambda0_1
  return( lambda0_1 = lambda0_1)
}

######################
#  estimating equation ########################
cox_equation_1 <- function(beta,matrix_proba_T,Ts,event,
                           X12,X34, p1,p2) {
  #X34=as.matrix(X34,ncol=p2,nrow=nB)
  # X12=as.matrix(X12,ncol=p1,nrow=nA)
  # p = length(p1)+length(p2)
  
 # eta1 <- X34 %*% beta[p2]        # length nB
  #eta2 <- X12 %*% beta[p1]    # length nA
  # form matrix eta_mat of size nA x nB: eta_mat[i,j] = eta1[i] + eta2[j]
#  eta_mat <- t(matrix(rep(as.numeric(eta1), times = nA), nrow = nB)) +
 #   t( matrix(rep(as.numeric(eta2), each  = nB), nrow = nB))
  
  eta_mat= expo_beta  (X12,X34, p1,p2,beta )
  
  ####### 
   pi_exp=  matrix_proba_T * exp(eta_mat)  # matrix nA*nB
  delta_ij <- matrix(event, nA, nB, byrow = TRUE)
 # W = matrix_proba_T*delta_ij
  
   RiskSets <- lapply(Ts, function(t) which(Ts >= t))
   
   H <- numeric(p)
   
   for (j in which(event == 1)) {
     
     rs <- RiskSets[[j]]
     denom <- sum(pi_exp[, rs])
     if (!is.finite( denom) || denom == 0)  denom <- .Machine$double.eps # avoid the 0
     
     #   numerator 
     EX34  <- colSums(pi_exp[, rs, drop=FALSE]) %*% X34[rs, , drop=FALSE]
     EX12 <- rowSums(pi_exp[, rs, drop=FALSE]) %*% X12
     
     EX <- c(EX12,EX34) / denom
     
     # equation 
     for (i in 1:nA) {
       Xij <- c(X12[i, ],X34[j, ])
       H <- H + matrix_proba_T[i, j] * (Xij - EX)
     }
   }
   return(H)

  #  H <- rep(0, p)
 
  #  for (i in 1:nA) {
   #   for (j in 1:nB) {
        
  #      X_ij = c (X1[j,], X12[i,])
 
  #risq = GetRiskSet (Ts[j], Ts)
 # nrisk = length(risq)
  
 # if(nrisk==1){
 #   denom_i = pi_exp_beta [,risq]  # denominator
 # }else {
 #   denom_i = rowSums( pi_exp_beta [,risq] )# sum on j
 # }
 #  denom  = sum (denom_i) #sum on i 
  # if (!is.finite( denom) || denom == 0)  denom <- .Machine$double.eps # avoid the 0
   
  # num_vec = rep(0, p)
   
  # for (k in 1:nA) {
  #   for (l in 1:nB) {
  #      if(Ts [j] <= Ts[l]){  
   #    x_kl <- c(X1[l, ], X12[k, ])
   #      num_vec <- num_vec + pi_exp_beta[k,l] * x_kl
  #     }
 #  }
  # }
  # H <- H + W [i, j] * (X_ij - num_vec / denom)
   #   }
    #  }
   # 
   # return(as.numeric(H))

}
   
#####################
# ###########solve the equation 

cox_1<- function( matrix_proba_T,Ts,event, X12,X34,
                  p1,p2,beta_ini,maxiter = 20){
  
  f <- function(x){
   val<-cox_equation_1( beta=x,matrix_proba_T,Ts,event, X12,X34,p1,p2)
   
   if (!is.numeric(val) || length(val) != length(x)) {
     stop(" cox_equation_1 must return a numeric vector of same length as x")
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

cox_itteration<-function(beta0,lambda0,Ts,event, X12,X34,
                         p1,p2, Q,tol= 1e-6, maxits = 300){
  
  X34=as.matrix(X34)
  X12=as.matrix(X12)
  nA=nrow(X12)
  nB=nrow(X34)
    p = length(p1)+length(p2)
 
  beta_ini=rep(0,p)
  it = 0
  converge = FALSE
  
  while ((!converge) && (it < maxits)){ 
    
    lambda0_old = lambda0
    beta0_old = beta0
    
    #expectation  
    
    lambda2 = Funct_lambda2(lambda0_old,Ts)
    
    matrix_proba_T = Func_prob_T (beta0_old ,lambda2, event, X12,X34,p1,p2, Q)
    
    
    #maximization
    lambda0 = Function_lambda0 (matrix_proba_T,beta0_old,Ts,event, X12,X34,p1,p2)
    
    estime = cox_1 ( matrix_proba_T,Ts,event, X12,X34, p1,p2,beta_ini,maxiter = 20)
    beta0 = estime$beta0
    beta_ini = beta0_old
    
    d_beta   <- sqrt(sum((beta0 - beta0_old)^2))
    d_lambda <- sqrt(sum((lambda0 - lambda0_old)^2))
    converge <- (d_beta + d_lambda) < tol
    #converge = sqrt (sum ((beta0.old -beta0))^2 + sum ( (lambda0.old - lambda0)^2 ) ) < tol  
    
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
  return(list(beta0=beta0,
              converge= converge, it=it))
}
