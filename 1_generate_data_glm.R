Generate_data <- function(beta,K,nA,nB,vector_pk,censor, min_prev = 0.01){
  
  # First database B with matching variables M
  datB = matrix(0, nrow = nB, ncol = K+1)
  
  conditionB = TRUE
  while (conditionB){
    datB[,1:K] = sapply(vector_pk, function(x){rbinom(nB, size = 1,prob = x)})
    conditionB = (sum(colSums(datB[,1:K]/nB) >= min_prev) < K)
  }
  datB[,K+1] = 1:nB #id
  
  datB = data.frame(datB)
  colnames(datB)=c(paste("M", 1:K, sep = ""),"id") 
  
  ##covarites data
  X1 = rnorm(nB,0,1)
  X2 = rbinom(nB,size = 1, prob = 0.7)
  X3 = rnorm(nB,0,2)
  X  = as.matrix(cbind(X1, X2, X3))
  
  #Survival data
  U  = runif(nB,0,1) # Time to event variable
  Tt = -log(U)/exp((X%*%beta))#lambda=1
  #c = quantile(Tt[,1], probs = 0.7)
  #c  = 1.777907
  c=  censor
  
  Time      = pmin(Tt,c) # le vrai temps
  delta     = as.numeric(Tt<=c)
  #X12 = cbind(X2,X3)
  #surv_data = data.frame(Time,delta,X1,X12)
  surv_data = data.frame(Time,delta,X1 ,X2,X3)
  
  #Y linear regression data
  varXbeta= beta[1]^2+ beta[2]^2*(0.7)*(0.3)+ beta[3]^2*4
   sigma1<- varXbeta * (1 - R2) / R2
  Y <- X%*%beta + rnorm(nB, mean = 0, sd = sqrt(sigma1))
 
  #Z classification data (logistic)
  proba_Z <- exp(X %*% beta) / (1 + exp(X %*% beta))
  US <- runif(nB)
  Z <- ifelse(US < proba_Z, 1, 0)# pourquoi pas rbin(nB,1,proba_Z)
  #Z<- rbinom(nB,1, proba_Z)
    
  #V poisson
  lambda_V <- exp(X %*% beta)
  V<- rpois(nB, lambda_V)
  
  ## matching variable A 
  idBA <- sample(1:nB,nA) #ident in A appearing in B
  datA = data.frame(datB)[idBA,]
    
## largest dabaset without match
  database = cbind.data.frame(Y, Z, V,Time,delta,X1,X2,X3,datB )
  
    return(list(database=database, datA=datA, datB=datB,
                 idBA=idBA) )
}
####################
#####################
#XB = cbind (X3)
#XA= cbind (X1,X2)
data_context <- function( K,idBA,Y, Z, V,Time,delta,XA,XB ,datB){
 
  # complete data (context 1) 
 databaseB_1 = cbind.data.frame(datB,Y, Z, V,Time,delta,XB ) [,-(K+1)]
databaseA_1 <- cbind.data.frame(datB,XA)[idBA,- (K+1) ]


# complete data (context 2)
databaseB_2 = cbind.data.frame(XB  )
databaseA_2 <- cbind.data.frame(Y, Z, V,Time,delta,XA)[idBA, ]

return( list( databaseB_1=databaseB_1,
              databaseA_1=databaseA_1,databaseB_2=databaseB_2,
              databaseA_2=databaseA_2))
}






