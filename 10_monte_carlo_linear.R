 
library(nleqslv)

setwd("C:/Users/fchezeut/Documents/GitHub/GLM/methods_glm")
source("1_generate_data_glm.R")
source("2_matrice_Q_inconnu.R")
source("2_1_function_exponential.R" )
source("3_naive_equation_glm.R")
source("6_Linear_glm.R")
source("8_Linear_method2.R")
source("3_scenarios_glm.R")


################### nsim monte carlos
estimates_parameters<- function(nsim,beta, K,nA,nB, V_A, V_B,
                                p1,p2,vector_pk,censor){
  p = length(p1)+length(p2)
  Linear_true_s = matrix(0,nrow = nsim, ncol = p) # continuous
  Linear_naive_s = matrix(0,nrow = nsim, ncol = p)
  
  coef_Linear_s = matrix(0,nrow = nsim, ncol = p)
  converge_Linear = vector()
  
  coef_Lin_method2_s = matrix(0,nrow = nsim, ncol = p)
  converge_Lin_method2 = vector()
  
  for (i in 1:nsim){
    
    #data generation
    data= Generate_data(beta,K,nA,nB,vector_pk,censor, min_prev = 0.01)
    
    datA=data$datA # matching variables
    datB=data$datB
    idBA=data$idBA
    database = data$database
    
    # # comparison process
    M= Matrix_Q ( K,datA,datB)
    Q=M$Q
    matrix_id=M$matrix_id
    
    # Data naive
    naive_id=M$naive_id
    
    Time = database$Time # survival variable
    delta=database$delta
    
    counting = database$Z # classification variable 
    
    continuous = database$Y# continuous variable 
    
    poisson = database$V # counting process
    
    # covariate
    
    XB =  cbind(database[, V_B ])   ## in B
    XA=  cbind(database[,  V_A ])   ## in A
    colnames(XB) = V_B
    colnames( XA) =  V_A
    
    context = data_context (K, idBA,Y=continuous, Z=counting , V=poisson,Time,delta,
                            XA,XB ,datB) 
    databaseA_1=context$databaseA_1 # base A with covariables
    databaseB_1=context$databaseB_1 # base B with outcomes
    ### true data
    true_data_glm = data.frame(Time=database$Time,delta=database$delta,X1_true=database$X1,
                               X2_true= database$X2,X3_true=database$X3,
                               Z=database$Z,Y=database$Y,V=database$V)[idBA,]
    
    data_naive=generate_naive_data (V_A, V_B, naive_id,databaseA_1, databaseB_1 )
    # explanatory naive variables
    X_1= as.matrix(cbind(data_naive$X1,data_naive$X2,data_naive$X3), nrow=nrow(data_naive), ncol= p  )
    Y_naive1=  as.matrix(data_naive$Y )
    
    X12 = cbind (databaseA_1[,V_A])
    X34= cbind (databaseB_1[,V_B])
    Y= databaseB_1$Y
    
    beta0=rep(0,p) # initialisation 
    mu2_0=rep(0, length(p2))
    beta_ini=rep(0,p )
    variance0 = 0.1
    variance1 = 0.1
    ### beta estimations
    beta_linear = linear_itteration(beta0,variance0, p1,p2, Y,
                                    X12, X34,Q, tol= 1e-6, maxits = 300)
    beta_linear_method2 =Linear_itt_method2 (beta0,mu2_0,variance0,variance1, p1,p2, Y, X12,X34, 
                                             matrix_id,   tol= 1e-6, maxits = 300) 
    
    f3= linear_naive_1(X_1,Y_naive1)
    
    lin = glm(Y ~ X1_true + X2_true  +X3_true ,
              family = gaussian,
              data = true_data_glm)
    ######################
    
    coef_Linear_s[i,]=as.vector(beta_linear$beta0)
    sigma_linear=beta_linear$sigma
    #conv_linear=beta_linear$converge
    it_linear=beta_linear$it
    
    coef_Lin_method2_s[i,]= as.vector(beta_linear_method2$beta)
    mu2 =as.vector(beta_linear_method2$mu2)
    #sigma_1=beta_linear_method2$variance1
    #sigma_0=beta_linear_method2$variance0
    #converge_Lin_method2 =beta_linear_method2$converge
    it_lin_2=beta_linear_method2$it 
    
    ########### naive estimation
    
    Linear_naive_s [i,]=as.vector(f3) 
    
    ################ true estimations
    
    Linear_true_s [i,] =as.vector(lin$coefficients)[-1]
    
  }
  
  return(list( Linear_true_s1=Linear_true_s[,1], Linear_true_s2=Linear_true_s[,2], Linear_true_s3=Linear_true_s[,3],
               Linear_naive_s1=Linear_naive_s[,1],Linear_naive_s2=Linear_naive_s[,2],Linear_naive_s3=Linear_naive_s[,3],
               coef_Linear_s1=coef_Linear_s[,1],coef_Linear_s2=coef_Linear_s[,2],coef_Linear_s3=coef_Linear_s[,3],
               coef_Lin_method2_s1= coef_Lin_method2_s[,1],coef_Lin_method2_s2= coef_Lin_method2_s[ ,2],coef_Lin_method2_s3= coef_Lin_method2_s[ ,3],           
               conv_linear=conv_linear,converge_Lin_method2=converge_Lin_method2 ) )
  } 
######################################

setwd("C:/Users/fchezeut/Documents/GitHub/GLM/Results_glm")

for (i in ( c(3,6,9) )){
 nsim=scenarios[i,1]
  K=scenarios[i,2]
  nA=scenarios[i,3]
  nB= scenarios[i,4]
  censor= scenarios[i,5]
  vector_pk = rep(prevalence_sample,K)
  
  results_sample = estimates_parameters (nsim,beta, K,nA,nB, V_A, V_B,
                                         p1,p2,vector_pk,censor) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/GLM/Results_linear_glm/","nsim=",nsim,
                    "_nA=",nA,"_nB=",nB, "_K=",K,"_censor=",censor,".Rdata")
  
  save(results_sample,file = filename)
}

