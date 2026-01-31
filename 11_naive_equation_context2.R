
generate_naive_data <- function(naive_id,databaseA_2, databaseB_2){

  #context 2
  
   X_naive2 <-  databaseB_2[naive_id,]$X3 
   colnames(X_naive2)= "X3"
   data_naive_2= cbind.data.frame(databaseA_2[,-c(1:K)],X_naive2)
   
  return( data_naive_2 = data_naive_2 )
} 

# equation naive COX context 2
# X_naive2 = as.matrix(data_naive_2[,-c(1:5)],nrow=nrow(data_naive_2),ncol=K)
# event2 = data_naive_2$delta
# Ts2 = data_naive_2$Time

 equa_naive_2 <- function(beta,X_naive2,event2,Ts2) {
  p=ncol(X_naive2)
 eXbeta = exp(X_naive2 %*% beta)
 XeXbeta = X_naive2*matrix(rep(exp(X_naive2 %*% beta),p), ncol = p) 
 n = length(Ts2)

 s = 0
 for(i in 1:n){ 
# at risk
  risq = which(( Ts2 >= Ts2[i]))
 if(length(risq)==1){

 num_naive = XeXbeta[risq,]
 denum_naive = eXbeta[risq] 

 }else if (length(risq)!=1){

   num_naive = colSums( XeXbeta[risq,])
  denum_naive = sum(eXbeta[risq])  }  

     s = s+ event2[i]* ( X_naive2[i,]- num_naive/denum_naive)
 }

  return(H_naive = s) ##naive estimator 
  }
##############
##############
# solve the naive equation context 2

 coxph_naive_2 <- function(beta_ini,Ts2,even2, X_naive2, maxiter = 20){

   p=ncol(X2)

 f <- function(x){
   equa_naive_2(beta =x,X_naive2,event2,Ts2)
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
 
 
 ########## naive logistic 
 #context 2
 
 #Z2 = data_naive_2$Z
 
 ##score_logistic_naive2 <- function(beta, X_naive2, Z2) {
 # Z_proba <- 1 / (1 + exp(-X_naive2 %*% beta))  # probability p_i
 #score <- t(X_naive2) %*% (Z2 - Z_proba)         # équation  score
 #return(score)
 #}
 ########## solve equation logistic
 
 #context 2
 #logistic_naive_2 <- function(beta_ini, X_naive2, Z2, maxiter = 20){
 
 # f <- function(x){
 #  score_logistic_naive2 (beta=x, X_naive2, Z2)
 #}
 #fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
 #beta0 <- as.numeric(fit_manual$x)
 #iterations <- fit_manual$iter
 #converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
 #if (iterations  == maxiter) {
 #  cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
 # converge = FALSE
 #}
 #return(list( beta0 = beta0, converge = converge, iterations = iterations) )
 #}  
 
 ######linear
 # context 2
 #Y2= as.matrix(data_naive_2$Y,ncol=p1, nrow=nA)
 
 #linear_naive_2<- function(X_naive2,Y2){
 # beta0= solve(t(X_naive2) %*% X_naive2, t(X_naive2) %*% Y2)
 # return(beta0)
 #}
 
 #############  naive Poisson estimating equations#########
 ##
 # context 2
 #V2 = as.matrix(data_naive_2$V, ncol=p1, nrow=nA)
 
 #score_poisson_naive2 <- function(beta, X_naive2, V2) {
 # lambda <- as.vector(exp(X_naive2 %*% beta))
 #score <- t(X_naive2) %*% (V2 - lambda)
 #as.vector(score)
 #return(score)
 #}
 
 # context 2
 #poisson_naive_2 <- function( beta_ini,X_naive2, V2, maxiter = 20){
 
 # f <- function(x){
 #  score_logistic_naive2 (beta=x, X_naive2, V2)
 #}
 #fit_manual <- nleqslv( beta_ini,f, method = "Newton" )
 #beta0 <- as.numeric(fit_manual$x)
 #iterations <- fit_manual$iter
 #converge <- as.numeric(( iterations< maxiter)& (fit_manual$scalex==1))
 #if (iterations  == maxiter) {
 # cat("WARNING! NOT CONVERGENT FOR NEWTON!", "\n")
 #converge = FALSE
 #}
 #return(list( beta0 = beta0, converge = converge, iterations = iterations) )
 #}  
 
 