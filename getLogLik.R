
#============================ define Loglikelihood starts
getLogLik <- function(x, E, comp=c(1:4), model=1){
  
  flags <- E$flags; lbs <- E$lbs; ubs <- E$ubs 
  J <- E$J; K <- E$K; A <- E$A; K1 <- E$K1; K0 <- E$K0; J2 <- E$J2; K2 <- E$K2; G <- E$G 

  # now translate parameter vector into individual parameters
  lambda_jk <- matrix(x[flags==1], J2, K2)
  
  q_jk <- exp(x[flags==2]) # CAUTION!!
  # 2016-8-25
  if(E$parametrizations[1]==1) q_jk <- matrix(q_jk,J,K1) #scalar
  if(E$parametrizations[1]==2) q_jk <- matrix(rep(q_jk,K1),J,K1)
  if(E$parametrizations[1]==3) q_jk <- matrix(rep(q_jk,each=J),J,K1)
  if(E$parametrizations[1]==4) q_jk <- matrix(q_jk, nrow=J) #matrix
  
  tau <- x[flags==3]
  eta <- exp(x[flags==4])
  
  Mtmp <- exp(x[flags==5])
  M_jak <- array(NA, c(K1,J,A)) #KJA to match Fs below
  if(E$parametrizations[2]==1) M_jak <- array(Mtmp, c(K1,J,A)) #KJAG
  if(E$parametrizations[2]==2) for(j in 1:J) M_jak[,j,] <- Mtmp[j]
  if(E$parametrizations[2]==3) for(a in 1:A) M_jak[,,a] <- Mtmp[a]
  if(E$parametrizations[2]==4) for(k in 1:K1) M_jak[k,,] <- Mtmp[k]
   
  # change 2016-10-6
  M_jak[,1,] <- M_jak[,1,]/12  
   
  d3 <- 1; if(E$parametrizations[3]==2) d3 <- A 
  H0 <- array(x[flags==6], c(K,K1-1,d3))
  R_kj <- matrix(exp(x[flags==7]), K, J)
  N0_ak <- matrix(exp(x[flags==8]), K, A-1)
  nu0 <- x[flags==9]
           
  Pi0 <- array(0, c(K,K1,d3))
  for(a in 1:d3){
    for(k in 1:K){
      Pi0[k,,a] <- exp(c(H0[k,,a],0))
      Pi0[k,,a] <- Pi0[k,,a]/sum(Pi0[k,,a])
    }
  }
  
  # intermediate parameteres 
  Pi <- array(0, c(K,K1,A,J))
  if(E$parametrizations[3]==1) for(a in 1:A) for(j in 1:J) Pi[,,a,j] <- Pi0[,,1] 
  if(E$parametrizations[3]==2) for(j in 1:J) Pi[,,,j] <- Pi0
  #Pi[,,a,j] <- rdirichlet(K, rep(1/K1, K1))
  #for(k in 1:K) {Pi[k,k,a,j] <- .9; Pi[k,-k,a,j]<-.1/(K-1)}
  
  fac <- 10^eta*exp(-tau*10)
  v_a <- E$agevec^eta*exp(-tau*E$agevec)/fac 
  Fs <- array(0, c(K1,J,A,G)) #KJAG
  for(i in 1:E$ncat) Fs[E$cat_KJG[i,1], E$cat_KJG[i,2], ,E$cat_KJG[i,3]] <- v_a * (q_jk[E$cat_KJG[i,2],E$cat_KJG[i,1]] * E$cat_effort[i]) #v_a*f_jgk  
  
  
  Fs_sumOverG <- apply(Fs, c(1:3),'sum')
  S <- exp(-(Fs_sumOverG + M_jak))
  
  
  #tmp <- matrix(NA,J,K); for(i in 1:nrow(dat$cat)) tmp[dat$cat$year[i],dat$cat$basin[i]] <- 1 #no NA
  
  
  nu <- array(nu0, c(G,J,K)) #this is CV
  
  
  loglik <- numeric(4)
  
 if(any(comp==1)){
  #============================ compute L1
  # next compute new columns for P for foo2 below 
  P <- E$rec_phi_ij*(E$rec_delta_gjk+(1-E$rec_delta_gjk)*lambda_jk[E$ind_rec_JK])
    
 # #-------------------- compute P using model 1    
  # add Delta, Pi and U part
  P <- P * Fs[E$ind_rec_KJAG] *(1-S[E$ind_rec_KJA])/(Fs_sumOverG[E$ind_rec_KJA]+M_jak[E$ind_rec_KJA])
  
  #P[which(lags==0)] <- P[which(lags==0)]*Pi[ind_rec_0]
  
  for(i0 in 1:E$nrec){
    i <- E$rec_IJKK1AG[i0,1]
    j <- E$rec_IJKK1AG[i0,2]
    k <- E$rec_IJKK1AG[i0,3]
    k1 <- E$rec_IJKK1AG[i0,4]
    a <- E$rec_IJKK1AG[i0,5]
    a1 <- min(a+E$rec_lags[i0], A)
    g <- E$rec_IJKK1AG[i0,6]
    
    #P[i0] <- P[i0]*( Fs[k1,j,a1,g]*(1-S[k1,j,a1])/(Fs_sumOverG[k1,j,a1]+M_jak) )  #U_a1gjk part
    
    if(E$rec_lags[i0]==0) P[i0] <- P[i0]*Pi[k,k1,a1,j] 
     
    if(E$rec_lags[i0]>0){
     for(j0 in 0:(E$rec_lags[i0]-1)){
       a0 <- min(a+j0,A)
       if(model==1){
         if(j0==0) Delta_iaj <- Pi[,,a0,i+j0] %*% diag(S[,i+j0,a0])
         else Delta_iaj <- Delta_iaj %*% (Pi[,,a0,i+j0] %*% diag(S[,i+j0,a0]) )
       }
       if(model==2){
         if(j0==0) Delta_iaj <- diag(as.vector(Pi[,,a0,i+j0] %*% S[,i+j0,a0]))
         else Delta_iaj <- diag(as.vector( Delta_iaj %*% (Pi[,,a0,i+j0] %*% S[,i+j0,a0]) ))       
       }
     }
     Delta_iaj <- Delta_iaj %*% Pi[,,min(a0+1, A), j] 
     P[i0] <- P[i0] * Delta_iaj[k,k1] 
    }
    
  }
  #--------------------------------------------
  #print(mean(P))  
     
#  par(mfrow=c(1,2), mar=c(2,2,.2,0)+.3,mgp=c(1.1,.2,0), tck=-0.01, cex.axis=.8, cex.lab=.9, cex.main=1)
#  hist(foo2$P1, xlab='Model 1',main='',prob=F, col='lightblue')
#  hist(foo2$P2, xlab='Model 2',main='',prob=F, col='lightgreen')
  
  
  #### # check if merging correctly
  #tmp <- foo2
  #tmp$gjk <-  paste(lev$rec[[4]][as.integer(substr(tmp$gjk,1,1))],lev$rec[[5]][as.integer(substr(tmp$gjk,3,3))], lev$rec[[6]][as.integer(substr(tmp$gjk,5,10L))] ,sep='_') 
  #write.csv(file='tmp.csv', tmp)
  
  tmp <- aggregate(P, list(E$rec_iak), 'sum'); names(tmp) <- c('iak','P')
  
  Rs <- cbind(E$rec_R,E$rec_delta_gjk,E$rec_phi_ij,lambda_jk[E$ind_rec_JK]); Ts <- E$tag_noRec; Ps <- P 
  ind_rec_KJAG <- E$ind_rec_KJAG
  save(file=paste('m',model,sep=''),Rs=Rs,Ts=Ts,Ps=Ps,ind_rec_KJAG=ind_rec_KJAG)  
    
  #foo3 <- merge(x=foo, y=foo3, by='iak', all=T)
  #foo3$P[is.na(foo3$P)] <- 0
  #foo3$RlogP[is.na(foo3$RlogP)] <- 0
  
  # avoid using merge
  #RlogP <- c(tmp$RlogP,0)[E$ind_rec2tag]
  tmpP <- c(tmp$P,0)[E$ind_rec2tag]
  
  #loglik[1] <- sum(RlogP + E$tag_noRec*log(1-P))
  loglik[1] <- sum(E$rec_R*log(P)) + sum(E$tag_noRec*log(1-tmpP))
  
 } 
  
  
 if(any(comp==2)){
  #============================ compute L2
  lambda <- lambda_jk[E$ind_cat_JK]
  delta <- E$cat_delta_gjk[E$ind_cat_positive]      
  
  loglik[2] <- sum( 
    E$cat_R_um[E$ind_cat_positive]*log(lambda) - E$cat_R_tot[E$ind_cat_positive]*log(delta + (1-delta)*lambda)
    )
 } 
  
  
 if(any(comp==3)){
  #============================ compute L_CT 
  N <- array(NA, c(K0,J,A))
  for(j in 1:J) for(k in 1:K0) N[k,j,1] <- R_kj[k,j]
  for(a in 2:A) for(k in 1:K0) N[k,1,a] <- N0_ak[k,a-1]
  
  for(a in 2:A)
   for(j in 2:J)
    for(k in 1:K0)
     N[k,j,a] <- sum(N[k,j-1,a-1]*Pi[k,,a-1,j-1]*S[,j-1,a-1])
  
  
  tmp <- E$cat_KJG
  EC <- array(0, c(K1,J,G,A))
  l_ct <- 0
  tmp$sEC <- rep(NA, nrow(tmp))
  for(i0 in 1:nrow(tmp)){
    j <- tmp$year[i0]
    k1 <- tmp$basin[i0]
    g <- tmp$gear[i0]
    
    EC[k1,j,g,] <- Fs[k1,j,,g]*(1-S[k1,j,])/(Fs_sumOverG[k1,j,]+M_jak[k1,j,]) * colSums(N[,j,]*Pi[1:K0,k1,,j])
    sEC <- log(sum(EC[k1,j,g,]))
    
    #VC <- (nu[g,j,k1]*sEC)^2
    
    # 2016-7-14
    VC <- log(nu[g,j,k1]^2 + 1)
    
    l_ct <- l_ct - 0.5*(log(VC) + (E$log_cat_catch[i0]-sEC)^2/VC)
    tmp$sEC[i0] <- sEC
  }
  
  loglik[3] <- l_ct
 } 
 
 
 
 
 if(any(comp==3)){ 
  #============================ compute L_CP
  
#  tmp <- dat$age
#  tmp$logE <- rep(NA, nrow(tmp))
#  for(i0 in 1:nrow(tmp)){
#   k1 <- tmp$rec_basin[i0]
#   j <- tmp$year[i0]
#   g <- tmp$gear[i0]
#   a <- tmp$agep[i0]
#   tmp$logE[i0] <- log( EC[k1,j,g,a]/sum(EC[k1,j,g,]) )
#   }
#   
#  loglik[4] <- sum( tmp$n_sum*tmp$p_perc*tmp$logE )  #l_cp
  
  sumEC <- apply(EC, c(1:3), 'sum')
  tmp <- EC[E$ind_age_KJGA]/sumEC[E$ind_age_KJG]
  inds <- which(tmp>0)
  
  loglik[4] <- sum( E$age_n_gjka[inds]*log(tmp[inds]) )  #l_cp 
  
  #sum( age_n_gjka*log(dat$age$p_perc) ) 
  
 } 
   
  return(loglik) 
  #return(-sum(loglik))
  
  #loglik <- as.data.frame(loglik)
  #row.names(loglik) <- c('l1_1','l1_2','l2') 
  #return(loglik)
}
#============================ define Loglikelihood starts