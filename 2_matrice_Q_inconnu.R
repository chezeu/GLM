
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
  
  matrix_id <- matrix(vect_gamma,nrow=nA,ncol = nB, byrow = TRUE)
 # # common variables
  number = apply(matrix_id,1, FUN= function(x){length(which(x==1))}) # number of individual
  
  Q = matrix_id/number
  
  naive_id = apply(Q,1, which.max)
  
  return(list(Q=Q, naive_id=naive_id))
  
}


library(matrixStats)

###############"" Record Linkage#######################

EM_binary <- function(  datA, datB, K, tol = 1e-5, maxits = 500){
 
   comp_mat = compare_binary(datA, datB, K)
   # Starting point
 
  nA <- nrow(datA)
  nB <-  nrow(datB)
  
  N <- nrow(comp_mat)
  p = nA/N
  
  prev = (colMeans(datB[,1:K]))
  
  u=rep(0.5,K)
  m=rep(0.5,K)
  #u = ((1-e)*(1-prev) + e*prev)*(1-prev) + prev*((1-e)*prev+e*(1-prev))
  #m = rep(1 - e,K)
  
  # initializations
  
  comp_mat  <- comp_mat[,1:K]
  
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converge = FALSE
  
  
  while ((!converge) & (it < maxits)){ 
    
    p.old = p
    m.old = m
    u.old = u
    ### E
    # Compute expectation
    
    m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
    u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
    
    
    probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
    probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
    
    
    g = p*probM/(p*probM+(1-p)*probU)
    
    ### Maximization
    g_mat = matrix(rep(g,K),ncol = K)
    
    p = sum(g)/N
    m = colSums(g_mat*comp_mat)/sum(g)
    u = colSums((1-g_mat)*comp_mat)/sum(1-g)
    
    if (length(which(m > 0.99999)) > 0) {
      m[which(m > 0.99999)] <- 0.99999
    }
    if (length(which(m < 1e-05)) > 0) {
      m[which(m < 1e-05)] <- 1e-05
    }
    if (length(which(u > 0.99999)) > 0) {
      u[which(u > 0.99999)] <- 0.99999
    }
    if (length(which(u < 1e-05)) > 0) {
      u[which(u < 1e-05)] <- 1e-05
    }
    
    
    it = it + 1
    
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) && 
      all(abs(u.old - u)/u.old < tol)
    
    if (it == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
  }
  if (length(which(g< 1e-05)) > 0) {
    g[which(g < 1e-05)] <- 0
  }
   g <- matrix(g,nrow=nA,ncol=nB, byrow = TRUE)
  
  return(list(g=g, p=p, m=m, u = u, it = it, converge = converge))
}

############## fellegi and sunter

FS_function = function(datA,datB,K){
# pairs with same blocking variable

datA$block = rep(1,nrow(datA))#variable suplementaire
datB$block = rep(1,nrow(datB))
paires <- pair_blocking(datA, datB,on="block",
                        add_xy = TRUE)
paires <- compare_pairs ( paires,
        on = paste0("M",1:K ),
         default_comparator = cmp_identical())
## classification
p0 = nA/ (nA*nB) # probability that a pair is a match
vars <- paste0("M", 1:K)

form <- as.formula(
  paste("~", paste(vars, collapse = " + "))
)

model <- problink_em( form, data = paires,  mprobs0 = list(0.9),
  uprobs0 = list(0.02), p0 = p0, tol = 1e-05,
  mprob_max = 0.999, uprob_min = 1e-04)

p_compar <- predict(model, paires, type ="mpost", add = TRUE, binary = TRUE)

prob_match <- matrix(p_compar$mpost, nrow = nA, ncol=nB,byrow = TRUE)

return( prob_match =prob_match )
}





