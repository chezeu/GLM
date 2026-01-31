
library(survival)
 
library(nleqslv)

setwd("C:/Users/fchezeut/Documents/GitHub/GLM/methods_glm")
source("1_generate_data_glm.R")
source("2_risk_function.R")
source("2_matrice_Q_inconnu.R")
source("3_naive_equation_glm.R")
source("4_Cox_glm.R")
source("5_Logistic_glm.R")
source("6_Linear_glm.R")
source("7_Poisson_glm.R")
source("8_Linear_method2.R")
source("3_scenarios_glm.R")

################### nsim monte carlos
estimates_parameters<- function(nsim, K,nA,nB, V_A,V_B,
                              p1,p2,vector_pk,censor){
  
  p=p1+p2
  
  Cox_true_s = matrix(0,nrow = nsim, ncol = p) # cox
  Cox_naive_s = matrix(0,nrow = nsim, ncol = p)
  coef_Cox_s = matrix(0,nrow = nsim, ncol = p)
  converge_Cox = vector()
  
  Logistic_true_s = matrix(0,nrow = nsim, ncol = p) # classification
  Logistic_naive_s = matrix(0,nrow = nsim, ncol = p)
  coef_Logistic_s = matrix(0,nrow = nsim, ncol = p)
  converge_Logistic = vector()
  
  Linear_true_s = matrix(0,nrow = nsim, ncol = p) # continuous
  Linear_naive_s = matrix(0,nrow = nsim, ncol = p)
  coef_Linear_s = matrix(0,nrow = nsim, ncol = p)
  converge_Linear = vector()
  
  Poisson_true_s = matrix(0,nrow = nsim, ncol = p) # counting
  Poisson_naive_s = matrix(0,nrow = nsim, ncol = p)
  coef_Poisson_s = matrix(0,nrow = nsim, ncol = p)
  converge_Poisson = vector()
  
  Lin_method2_true_s = matrix(0,nrow = nsim, ncol = p) # counting
  Lin_method2_naive_s = matrix(0,nrow = nsim, ncol = p)
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
   
   data_naive=generate_naive_data (V_A, V_B, naive_id,databaseA_1, databaseB_1 )
 
    Ts=databaseB_1$Time
    event=databaseB_1$delta
     X12 = cbind (databaseA_1[,V_A])
     X34= cbind (databaseB_1[,V_B])
     Z= databaseB_1$Z
     Y= databaseB_1$Y
     V=databaseB_1$V
     
     
    beta0=rep(0,p1+p2) # initialisation 
    mu2_0=rep(0, p2 )
    lambda0 =  rep(0.1,length(event))
    lambda0 [which(event == 0)] = 0
    beta_ini=rep(0,p1+p2)
    variance0 = 0.1
    variance1 = 0.1
    ### beta estimations
    beta_cox = cox_itteration (beta0,lambda0,Ts,event, X12,X34,p1,p2, Q,tol= 1e-6,
                               maxits = 300)# too slow
    beta_logistic = logistic_itteration(beta0, p1,p2, Z,
                                         X12,X34, Q, tol= 1e-6, maxits = 300)
    beta_linear = linear_itteration(beta0,variance0, p1,p2, Y,
                                   X12, X34,Q, tol= 1e-6, maxits = 300)
    beta_Poisson = Poisson_itteration (beta0,V, X12,X34, p1,p2, Q,tol= 1e-6,
                                   maxits = 300) #too slow
    beta_linear_method2 =Linear_itt_method2 (beta0,mu2_0,variance0,variance1, p1,p2, Y, X12,X34, 
                                             matrix_id,   tol= 1e-6, maxits = 300) 
      
    
    coef_Cox_s[i,] = as.vector(beta_cox$beta0)
    conv_cox=beta_cox$converge
    it_cox=beta_cox$it
    
    coef_Logistic_s[i,]= as.vector(beta_logistic$beta0)
    conv_logistic=beta_logistic$converge
    it_logistic=beta_logistic$it
    
    coef_Linear_s[i,]=as.vector(beta_linear$beta0)
    sigma_linear=beta_linear$sigma
    conv_linear=beta_linear$converge
    it_linear=beta_linear$it
   
    coef_Poisson_s[i,]= as.vector(beta_Poisson$beta0)
    conv_Poisson=beta_Poisson$converge
    it_Poisson=beta_Poisson$it
    
    coef_Lin_method2_s[i,]= as.vector(beta_linear_method2$beta)
    converge_Lin_method2 =beta_linear_method2$converge
    it_lin_2=beta_linear_method2$it 
    
    ########### naive estimation
    databaseA_2=data$databaseA_2 # CONTEXT 2 to be remove
     databaseB_2=data$databaseB_2
     # naive data
    data_naive_1=data_naive$data_naive_1
    
    T_naive1=   data_naive_1$Time 
    delta_naive =  data_naive_1$delta
    Z_naive1=  as.matrix(data_naive_1$Z,ncol=p1, nrow=nA)
    Y_naive1=  as.matrix(data_naive_1$Y,ncol=p1, nrow=nA)
    V_naive1=  as.matrix(data_naive_1$V,ncol=p1, nrow=nA)
    
   # explanatory variables
    X_1= as.matrix(cbind(data_naive_1$X1,data_naive_1$X2,data_naive_1$X3), nrow=nrow(data_naive_1), ncol= p  )
    
    #####
    f1= coxph_naive_1(beta_ini,Ts1=T_naive1,event1=delta_naive, X_1, maxiter = 20)
    Cox_naive_s[i,]= as.vector(f1$beta0)
    
   f2= logistic_naive_1 ( beta_ini=beta0, X_1, Z_naive1, maxiter = 20)
   Logistic_naive_s [i,]=as.vector(f2$beta0) 
    
    f3= linear_naive_1(X_1,Y_naive1)
    Linear_naive_s [i,]=as.vector(f3) 
    
    f4 = poisson_naive_1 (beta_ini, X_1, V_naive1, maxiter = 20)
    Poisson_naive_s[i,]= as.vector(f4$beta0) 
    
    ################ true estimations
  
   
    ### true data cox
    true_data_cox = data.frame(Time_true=database$Time,delta_true=database$delta,X1_true=database$X1,
                               X2_true= database$X2,X3_true=database$X3 )[idBA,]

    fit_true = coxph(Surv(Time_true,delta_true)~.,data = true_data_cox)
    Cox_true_s [i,]= as.vector(fit_true$coefficients)
    
    ### true logistic
    true_data_glm = data.frame(Time=database$Time,delta=database$delta,X1_true=database$X1,
                   X2_true= database$X2,X3_true=database$X3,Z=database$Z,Y=database$Y,V=database$V)[idBA,]
     
    logist = glm(Z ~ X1_true + X2_true  +X3_true ,
                 family = binomial(link = "logit"),
                 data = true_data_glm)
    Logistic_true_s [i,] = as.vector(logist$coefficients)[-1]
 
    ### true linear
    lin = glm(Y ~ X1_true + X2_true  +X3_true ,
                 family = gaussian,
                 data = true_data_glm)
    Linear_true_s [i,] =as.vector(lin$coefficients)[-1]

    ### true poisson
    Poi = glm(V ~ X1_true + X2_true  +X3_true ,
                   family = poisson(link = "log"),
                   data = true_data_glm)
    Poisson_true_s [i,] = as.vector(Poi$coefficients)[-1]
  
  }
  
  return(list( Cox_true_s1 = Cox_true_s[,1],Cox_true_s2 = Cox_true_s[,2],Cox_true_s3 = Cox_true_s[,3],
       Logistic_true_s1 = Logistic_true_s[,1],Logistic_true_s2 = Logistic_true_s[,2], Logistic_true_s3 = Logistic_true_s[,3],
       Linear_true_s1=Linear_true_s[,1], Linear_true_s2=Linear_true_s[,2], Linear_true_s3=Linear_true_s[,3],
       Poisson_true_s1 = Poisson_true_s[,1], Poisson_true_s2 = Poisson_true_s[,2], Poisson_true_s3 = Poisson_true_s[,3],
       Cox_naive_s1=Cox_naive_s[,1],Cox_naive_s2=Cox_naive_s[,2],Cox_naive_s3=Cox_naive_s[,3],
       Logistic_naive_s1=Logistic_naive_s[,1],Logistic_naive_s2=Logistic_naive_s[,2],Logistic_naive_s3=Logistic_naive_s[,3],
      Linear_naive_s1=Linear_naive_s[,1],Linear_naive_s2=Linear_naive_s[,2],Linear_naive_s3=Linear_naive_s[,3],
      Poisson_naive_s1=Poisson_naive_s[,1],Poisson_naive_s2=Poisson_naive_s[,2],Poisson_naive_s3=Poisson_naive_s[,3],
       coef_Cox_s1=coef_Cox_s[,1],coef_Cox_s2=coef_Cox_s[,2], coef_Cox_s3=coef_Cox_s[,3],
      coef_Logistic_s1=coef_Logistic_s[,1],coef_Logistic_s2=coef_Logistic_s[,2],coef_Logistic_s3=coef_Logistic_s[,3],
      coef_Linear_s1=coef_Linear_s[,1],coef_Linear_s2=coef_Linear_s[,2],coef_Linear_s3=coef_Linear_s[,3],
      coef_Poisson_s1=coef_Poisson_s[,1],coef_Poisson_s2=coef_Poisson_s[,2],coef_Poisson_s3=coef_Poisson_s[,3],
      conv_cox= conv_cox,conv_logistic=conv_logistic,
        conv_linear=conv_linear,conv_Poisson=conv_Poisson) ) } 
######################################

setwd("C:/Users/fchezeut/Documents/GitHub/GLM/Results_glm")

for (i in (1:nrow(scenarios))){
  nsim=scenarios[i,1]
  K=scenarios[i,2]
  nA=scenarios[i,3]
  nB= scenarios[i,4]
  censor= scenarios[i,5]
  vector_pk = rep(prevalence_sample,K)
  
  results_sample = estimates_parameters (nsim,beta0,K,nA,nB, 
                                         p1,p2,vector_pk,censor) #monte carlos
  
  filename = paste0("C:/Users/fchezeut/Documents/GitHub/GLM/Results_glm/","nsim=",nsim,
                    "_nA=",nA,"_nB=",nB, "_K=",K,"_censor=",censor,".Rdata")
  
  save(results_sample,file = filename)
}

