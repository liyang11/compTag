
translate <- function(x,E,model=1){
  flags <- E$flags; lbs <- E$lbs; ubs <- E$ubs 
  J <- E$J; K <- E$K; A <- E$A; K1 <- E$K1; K0 <- E$K0; J2 <- E$J2; K2 <- E$K2; G <- E$G 

  lambda_jk <- matrix(x[flags==1], J2, K2)
  out <- list(lambda_jk)
  names(out) <- 'lambda_jk'
  out$q_jk <- exp(x[flags==2])
  out$tau <- x[flags==3]
  out$eta <- exp(x[flags==4])
  out$M_jak <- exp(x[flags==5])
  
  H0 <- matrix(x[flags==6], K, K1-1)
  Pi0 <- matrix(0, K, K1)
  for(k in 1:K){
    Pi0[k,] <- exp(c(H0[k,],0))
    Pi0[k,] <- Pi0[k,]/sum(Pi0[k,])
  }
  out$Pi <- Pi0
  
  out$R_kj <- matrix(exp(x[flags==7]), K, J)
  out$N0_ak <- matrix(exp(x[flags==8]), K, A-1)
  
  out$nu0 <- x[flags==9]
  
  q_jk <- matrix(out$q_jk, J,K1)
  
  fac <- 10^out$eta*exp(-out$tau*10)
  v_a <- E$agevec^out$eta*exp(-out$tau*E$agevec)/fac 
  Fs <- array(0, c(K1,J,A,G)) #KJAG
  for(i in 1:E$ncat) Fs[E$cat_KJG[i,1], E$cat_KJG[i,2], ,E$cat_KJG[i,3]] <- v_a * (q_jk[E$cat_KJG[i,2],E$cat_KJG[i,1]] * E$cat_effort[i]) #v_a*f_jgk  
  
  
  Fs_sumOverG <- apply(Fs, c(1:3),'sum')
  S <- exp(-(Fs_sumOverG + out$M_jak))
  
  out$Fs <- Fs
  out$Fs_sumOverG <- Fs_sumOverG
  out$S <- S
  
  source('getLogLik.R')
  out$logLik <- getLogLik(x,E,comp=1:4,model=model)
  
  return(out)
}