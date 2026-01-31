
expo_beta <- function(X12,X34, p1,p2 ){
  
  #X3=as.matrix(X3,ncol=p2,nrow=nB)
  # X12=as.matrix(X12,ncol=p1,nrow=nA),X12=cbind(databaseA_1$X1,databaseA_1$X1)
  #nA=nrow(X12)
  # nB=nrow(X3)
  # p = p1+p2
  
  eta1 <- X34 %*% beta[-c(1:p1)]        # length nB
  eta2 <- X12 %*% beta[1:p1]    # length nA
  # form matrix eta_mat of size nA x nB: eta_mat[i,j] = eta1[i] + eta2[j]
  eta_mat <- t(matrix(rep(as.numeric(eta1), times = nA), nrow = nB)) +
    t( matrix(rep(as.numeric(eta2), each  = nB), nrow = nB))
  
  return(eta_mat=eta_mat)
}

