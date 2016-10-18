
#============================ define Loglikelihood starts
getLogLik <- function(x, E, comp=c(1:4), model=1){
  
  flags <- E$flags; lbs <- E$lbs; ubs <- E$ubs 
  I <- E$I; J <- E$J; K <- E$K; A <- E$A; K1 <- E$K1; K0 <- E$K0; J2 <- E$J2; K2 <- E$K2; G <- E$G

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
  nu <- array(nu0, c(G,J,K)) #this is CV
  
  loglik <- numeric(4)
  
  if(any(comp==1)){ #************  compute L1 starts

   # construct P-array
   P <- array(0, c(I,A,K,G,J,K1))
   phiFun <- function(x) return(1-0.31/(1+exp(8.52-x)))
   
   #******** We don't need to calculate P at iak where there is no released fish in real data because R and T_noRecapture will be both 0, no contribution to likelihood
   # cat_gjk_raw <- paste(E$cat_KJG$gear,E$cat_KJG$year,E$cat_KJG$basin,sep='_')
#   for(i in 1:I){for(a in 1:A){for(k in 1:K){ #loop iak starts
#    ind1 <- which(E$tag_iak==paste(i,a,k,sep='_'))
#    if(length(ind1)>0){ 
    
   for(m0 in 1:E$ntag){ # this is equivalent to for(i in 1:I){for(a in 1:A){for(k in 1:K) above 
     i <- E$tag_dbIAK[m0,1]
     a <- E$tag_dbIAK[m0,2]
     k <- E$tag_dbIAK[m0,3] 
     
     #******** We don't need to calculate P at gjk where there is no fishing efforts in real data, because catch_effort=0 -> U=0 -> P=0
#     for(g in 1:G){for(j in i:J){for(k1 in 1:K1){#loop gjk starts
#      ind2 <- which(cat_gjk_raw==paste(g,j,k1,sep='_')) 
#      if(length(ind2)>0){ 
    for(m1 in 1:E$ncat){ # this is equivalent to for(g in 1:G){for(k1 in 1:K1){for(j in i:J) above
      g <- E$cat_KJG[m1,3] 
      j <- E$cat_KJG[m1,2]  
      k1 <- E$cat_KJG[m1,1]
      if(j>=i){
      
       a1 <- min( a+ (j-i), A )
       i1 <- 1 #if(j<=4) i1 <- 1 else i1 <- 2
       if(k1<=3) i2 <- 1 else if(k1<=5) i2 <- 2 else if(k1<=8) i2 <- 3 else if(k1==9) i2 <- 4 else i2 <- 5
       delta_gjk0 <- E$cat_delta_gjk[m1] #[ind2]
       phi_ij <- phiFun( (j - i)*12 + 1 )
       P[i,a,k,g,j,k1] <- phi_ij*(delta_gjk0+(1-delta_gjk0)*lambda_jk[i1,i2])
       P[i,a,k,g,j,k1] <- P[i,a,k,g,j,k1] * Fs[k1,j,a1,g] *(1-S[k1,j,a1])/(Fs_sumOverG[k1,j,a1]+M_jak[k1,j,a1])  
         
     if(j==i) P[i,a,k,g,j,k1] <- P[i,a,k,g,j,k1]*Pi[k,k1,a1,j] 
     if(j>i){
         lags0 <- j-i
     for(j0 in 0:(lags0-1)){
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
     P[i,a,k,g,j,k1] <- P[i,a,k,g,j,k1] * Delta_iaj[k,k1] 
    }
             
   #} }}} ##loop gjk ends
   }}
   
 #} }}}  #loop iak ends 
 }
 
  tmpP = apply(P, c(1:3),'sum')
  loglik[1] <- sum(E$rec_R*log(P[E$ind_rec_i6])) + sum(E$tag_noRec*log(1-tmpP[E$ind_tag_i6]))
  
 } #************  compute L1 ends
 
  
 if(any(comp==2)){ #************  compute L2 starts
  lambda <- lambda_jk[E$ind_cat_JK]
  delta <- E$cat_delta_gjk[E$ind_cat_positive]     
  loglik[2] <- sum(  E$cat_R_um[E$ind_cat_positive]*log(lambda) - E$cat_R_tot[E$ind_cat_positive]*log(delta + (1-delta)*lambda) )
 } #************  compute L2 ends
  
  
 if(any(comp==3)){ #************  compute L_CT starts
  
  N <- array(NA, c(K0,J,A))
  for(j in 1:J) for(k in 1:K0) N[k,j,1] <- R_kj[k,j]
  for(a in 2:A) for(k in 1:K0) N[k,1,a] <- N0_ak[k,a-1]
  
  for(a in 2:A)
   for(j in 2:J)
    for(k in 1:K0)
     N[k,j,a] <- sum(N[k,j-1,a-1]*Pi[k,,a-1,j-1]*S[,j-1,a-1])
  
  EC <- array(0, c(K1,J,G,A))
  tmp <- E$cat_KJG
  l_ct <- 0
  for(i0 in 1:nrow(tmp)){
    j <- tmp$year[i0]
    k1 <- tmp$basin[i0]
    g <- tmp$gear[i0]
    EC[k1,j,g,] <- Fs[k1,j,,g]*(1-S[k1,j,])/(Fs_sumOverG[k1,j,]+M_jak[k1,j,]) * colSums(N[,j,]*Pi[1:K0,k1,,j])
    sEC <- log(sum(EC[k1,j,g,]))
    VC <- log(nu[g,j,k1]^2 + 1)
    l_ct <- l_ct - 0.5*(log(VC) + (E$log_cat_catch[i0]-sEC)^2/VC)
  }
  loglik[3] <- l_ct
  
 }#************  compute L_CT ends 
 
 
 if(any(comp==3)){ #************  compute L_CP starts
  sumEC <- apply(EC, c(1:3), 'sum')
  tmp <- EC[E$ind_age_KJGA]/sumEC[E$ind_age_KJG]
  inds <- which(tmp>0)
  loglik[4] <- sum( E$age_n_gjka[inds]*log(tmp[inds]) )
 } #************  compute L_CP ends
   
 return(loglik)  
}
#============================ define Loglikelihood ends