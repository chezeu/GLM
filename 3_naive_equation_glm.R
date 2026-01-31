 
generate_naive_data <- function(V_A, V_B, naive_id,databaseA_1, databaseB_1){
                                
  #context 1
 X12_naive = cbind (databaseA_1[,V_A])
 X34_naive = cbind( databaseB_1[naive_id,V_B])
 colnames(X12_naive)= V_A
 colnames(X34_naive)= V_B
  Y_naive <- databaseB_1[naive_id ,c("Y","Z","V","Time","delta")]
  data_naive_1= cbind.data.frame( Y_naive, X12_naive,
                                  X34_naive)
 
  return( data_naive_1 = data_naive_1 )
}

#################
# equation naive COX context 1

equa_naive_1 <- function(beta,X_1,event1,Ts1) {
  p=ncol(X_1)
  eXbeta = exp(X_1 %*% beta)
  XeXbeta = X_1*matrix(rep(exp(X_1 %*% beta),p), ncol = p) 
  n = length(Ts1)
  
  s = 0
  for(i in 1:n){ 
    # at risk
    risq = which(( Ts1 >= Ts1[i]))
    if(length(risq)==1){
      
      num_naive = XeXbeta[risq,]
      denum_naive = eXbeta[risq] 
      
    }else if (length(risq)!=1){
      
      num_naive = colSums( XeXbeta[risq,])
      denum_naive = sum(eXbeta[risq])  }  
    
    s = s+ event1[i]* ( X_1[i,]- num_naive/denum_naive)
  }
  # s
  return(H_naive = s) ##naive estimator 
}

# solve the naive equation context 1

coxph_naive_1 <- function(beta_ini,Ts1,event1, X_1, maxiter = 20){
  
  p=ncol(X_1)
  
  f <- function(x){
    equa_naive_1(beta =x,X_1,event1,Ts1)
  }
  fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list(beta0 = beta0, converge = converge, iterations = iterations))
}

 ######### Logistic ########################################
###########################################
# naive loistic equation
#context 1
#Z_naive1=  data_naive_1$Z

score_logistic_naive1 <- function(beta, X_1, Z_naive1) {
  Z_proba <- 1 / (1 + exp(-X_1 %*% beta))  # probability p_i
  score <- t(X_1) %*% (Z_naive1 - Z_proba)         # équation  score
  return(score)
}

########## solve equation logistic
#context 1
logistic_naive_1 <- function( beta_ini, X_1, Z_naive1, maxiter = 20){
  
  f <- function(x){
    score_logistic_naive1(beta=x, X_1, Z_naive1)
  }
  fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
   if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}
####### naive linear estimation#################
##############################################################33

# context 1
#Y_naive1 = as.matrix(data_naive_1$Y,ncol=p1, nrow=nA)

linear_naive_1<- function(X_1,Y_naive1){
  beta0 = solve(t(X_1) %*% X_1, t(X_1) %*% Y_naive1)
  return(beta0)
}

###############  naive Poisson estimating equations#########
###############################################################
# context 1
#V_naive1 = as.matrix(data_naive_1$V, ncol=p1, nrow=nA)

score_poisson_naive1 <- function(beta, X_1, V_naive1) {
  lambda <- as.vector(exp(X_1 %*% beta))
  score <- t(X_1) %*% (V_naive1 - lambda)
  as.vector(score)
  return(score)
}

################## solve the poisson estimating equations
# context 1
poisson_naive_1 <- function(beta_ini, X_1, V_naive1, maxiter = 20){
  
  f <- function(x){
    score_logistic_naive1 (beta=x, X_1, V_naive1)
  }
  fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
  beta0 <- as.numeric(fit_manual$x)
  iterations <- fit_manual$iter
  converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
  if (iterations  == maxiter) {
    cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
    converge = FALSE
  }
  return(list( beta0 = beta0, converge = converge, iterations = iterations) )
}  



