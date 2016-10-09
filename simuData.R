

#============================ simulate data using para
# now not fix p(observed R =0) = 0
simuData <- function(x, E, seeds){
  set.seed(seeds)
  
  model <- 2 
 
  flags <- E$flags; lbs <- E$lbs; ubs <- E$ubs
  I <- E$I; G <- E$G 
  J <- E$J; K <- E$K; A <- E$A; K1 <- E$K1; K0 <- E$K0; J2 <- E$J2; K2 <- E$K2; G <- E$G 
  
  loglik <- numeric(4) # store true likelihood
  
  # now translate parameter vector into individual parameters
  # now translate parameter vector into individual parameters
  lambda_jk <- matrix(x[flags==1], J2, K2)
  
  q_jk <- exp(x[flags==2]) # CAUTION!!
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
  for(i in 1:E$ncat) Fs[E$cat_KJG[i,1], E$cat_KJG[i,2], ,E$cat_KJG[i,3]] <- v_a * (q_jk[E$cat_KJG[i,2],E$cat_KJG[i,1]] * E$cat_effort[i]) #note we are using observed cat_effort!!
  
  
  Fs_sumOverG <- apply(Fs, c(1:3),'sum')
  S <- exp(-(Fs_sumOverG + M_jak))

  nu <- array(nu0, c(G,J,K)) #this is CV
  
  phiFun <- function(x) return(1-0.31/(1+exp(8.52-x)))
  phi_ij <- matrix(NA,I,J)
  for(i in 1:I)
   for(j in 1:J)
    phi_ij[i,j] <- phiFun( (j - i)*12 + 1 )
  
  
  N <- array(NA, c(K0,J,A))
  for(j in 1:J) for(k in 1:K0) N[k,j,1] <- R_kj[k,j]
  for(a in 2:A) for(k in 1:K0) N[k,1,a] <- N0_ak[k,a-1]
  
  for(a in 2:A)
   for(j in 2:J)
    for(k in 1:K0)
     N[k,j,a] <- sum(N[k,j-1,a-1]*Pi[k,,a-1,j-1]*S[,j-1,a-1])
  
  
  catch_m <- exp(E$log_cat_catch)*E$cat_delta_gjk # note this is the actually observed catch_m before we overwrite E$log_cat_catch below
  #tmp <- dat$cat
  EC <- array(0, c(K1,J,G,A))
  l_ct <- 0
#  cat_catch <- array(0, c(J,K1,G))
  log_cat_catch <- numeric(nrow(E$cat_KJG))
  #age_n_gjka <- array(0, c(J,K1,G,A)) #n_gjk'a
  age_n_gjka <- ind_age_KJG <- ind_age_KJGA <- numeric()
#  tmp <- numeric() 
#  for(j in 1:J){
#   for(k1 in 1:K1){
#    for(g in 1:G){      
#     if(any(Fs[k1,j,,g]!=0)){ 
   # there's no need to go through this 4 conditions since we are using the real catch efforts. It would match with E$cat_KJG. So this would be the same as likelihood calculation based on real catch efforts stored in E$cat_KJG, therefore I just copied:
  tmp <- E$cat_KJG
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
    #************** simulate log-catch from normal truncated below
    #log_cat_catch[i0] <- catch_m[i0]-1; while(log_cat_catch[i0]<catch_m[i0]) log_cat_catch[i0] <- exp(rnorm(1,sEC,sqrt(VC)))
    log_cat_catch[i0] <- rnorm(1,sEC,sqrt(VC))
    if(log_cat_catch[i0]<log(catch_m[i0])) log_cat_catch[i0] <- rtruncnorm(1,sEC,sqrt(VC),a=log(catch_m[i0]),b=Inf)
    l_ct <- l_ct - 0.5*(log(VC) + (log_cat_catch[i0]-sEC)^2/VC)
    
    #**************  simulate age composition age_n_gjka
    p_perc <- EC[k1,j,g,]/sum(EC[k1,j,g,])
    #simulate age_n_gjka
    #age_n_gjka[j,k1,g,] <- 25*p_perc # not generated from MN?
    vec <- rmultinom(1,25, prob=p_perc)
    age_n_gjka <- c(age_n_gjka, vec)
    
    # vectorization
    foo <- array(0, c(K,J,G)); foo[k1,j,g] <- 1
    ind_age_KJG <- c(ind_age_KJG, rep(which(foo==1), A))
    for(a in 1:A){
      foo <- array(0, c(K,J,G,A)); foo[k1,j,g,a] <- 1
      ind_age_KJGA <- c(ind_age_KJGA,which(foo==1))
    }
  }
#      tmp <- rbind(tmp, c('k1','j','g',log(cat_catch[j,k1,g]),cat_delta_gjk))
#     }
#    }
#   }
#  }
  loglik[3] <- l_ct
  
  Cs <- cbind(E$cat_KJG, catch_m, E$cat_effort, exp(E$log_cat_catch),E$cat_delta_gjk, E$cat_R_tot, E$cat_R_um)
  names(Cs) <- c('k1','j','g','catch_m','efforts','cat_catch_real','delta_gjk_real','R_tot_real','R_um_real')
  
  
  #********************  overwrite session 1 ******************
  #cbind(catch_m, exp(log_cat_catch))
  E$cat_delta_gjk <- catch_m/exp(log_cat_catch) # note we are using observed catch_m
  which(E$cat_delta_gjk>1)
  if(any(E$cat_delta_gjk>1)) stop('simulated catch is less than the observed catch_m')
  E$log_cat_catch <- log_cat_catch
  E$age_n_gjka <- age_n_gjka
  E$ind_age_KJG <- ind_age_KJG
  E$ind_age_KJGA <- ind_age_KJGA
  #************************************************************
  
  # we are not finishing updating simulated data yet as we need also simulate R_um from Likelihood 2: Binomial(R_tot, p), but before that, we will need first simulate R_tot from Likelihood 1.  
  
  Cs$cat_catch_simu <- exp(E$log_cat_catch)
  Cs$delta_gjk_simu <- E$cat_delta_gjk
  
  sumEC <- apply(EC, c(1:3), 'sum')
  loglik[4] <- sum( age_n_gjka*log(EC[ind_age_KJGA]/sumEC[ind_age_KJG]) )  #l_cp 
   
  
  #============================ compute L1
  # next compute new columns for P for rec_cat below 
  P <- array(0, c(I,A,K,G,J,K1))
  Rs <- Ts <- Ps <- numeric()
  cat_gjk_raw <- paste(E$cat_KJG$gear,E$cat_KJG$year,E$cat_KJG$basin,sep='_')
  
  
  kk <- 1; ind_rec_KJAG <- integer()
  
  for(i in 1:I){for(a in 1:A){for(k in 1:K){ #loop iak starts
    ind1 <- which(E$tag_iak==paste(i,a,k,sep='_'))
    if(length(ind1)>0){ # We don't simulate at iak where there is no released fish in real data
     delta_gjkTmp <- array(NA,c(G,J,K1))
     for(g in 1:G){for(j in i:J){for(k1 in 1:K1){#loop gjk starts
      a1 <- min( a+ (j-i), A )
      ind2 <- which(cat_gjk_raw==paste(g,j,k1,sep='_')) 
      if(length(ind2)>0){ # We don't simulate at gjk where there is no fishing efforts in real data
       i1 <- 1 #if(j<=4) i1 <- 1 else i1 <- 2
if(k1<=3) i2 <- 1 else if(k1<=5) i2 <- 2 else if(k1<=8) i2 <- 3 else if(k1==9) i2 <- 4 else i2 <- 5
       #delta_gjk0 <- mean(dat$cat$catch_m/dat$cat$catch) 
       #use observed effort, catch_m from our true data for the simulation study.
       #delta_gjk0 <- 0; inds <- which(dat$cat$gear==g&dat$cat$year==j&dat$cat$basin==k1); if(length(inds)>0){if(length(inds)==1) {delta_gjk0 <- dat$cat$catch_m[inds]/dat$cat$catch[inds]} else stop('multiple match')} 
       delta_gjkTmp[g,j,k1] <- delta_gjk0 <- E$cat_delta_gjk[ind2]
       P[i,a,k,g,j,k1] <- phi_ij[i,j]*(delta_gjk0+(1-delta_gjk0)*lambda_jk[i1,i2])
       P[i,a,k,g,j,k1] <- P[i,a,k,g,j,k1] * Fs[k1,j,a1,g] *(1-S[k1,j,a1])/(Fs_sumOverG[k1,j,a1]+M_jak[k1,j,a1]) 
         
     a1 <- min( a+ j-i, A )
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
         
   } }}} ##loop gjk ends
   # simulate using observed total release tag_rec$T
   nn <- length(tmp <- rmultinom(1,E$tag_tot[ind1],prob=c(as.numeric(P[i,a,k,,,]),1-sum(P[i,a,k,,,])) ))
   Rtmp <-  array(tmp[1:(nn-1)], c(G,J,K1)) #first nn-1 are the simulated R(recapture), the last is tag_noRec = tot-Riak...
   for(g in 1:G){
      for(j in 1:J){
       for(k1 in 1:K1){
        if(Rtmp[g,j,k1]>0){
          i1 <- 1 #if(j<=4) i1 <- 1 else i1 <- 2
if(k1<=3) i2 <- 1 else if(k1<=5) i2 <- 2 else if(k1<=8) i2 <- 3 else if(k1==9) i2 <- 4 else i2 <- 5
          Rs <- rbind(Rs, c(i,a,k,g,j,k1,Rtmp[g,j,k1],phiFun((j-i)*12+1),delta_gjkTmp[g,j,k1],lambda_jk[i1,i2]))
          
          a1 <- min( a+ (j-i), A )
          tmp1 <- array(0, c(K1,J,A,G)); tmp1[k1,j,a1,g] <- 1
          ind_rec_KJAG[kk] <- which(tmp1==1); kk <- kk+1
          
          Ps <- c(Ps, P[i,a,k,g,j,k1])
          }
       }}}
 Ts <- rbind(Ts, c(i,a,k,E$tag_tot[ind1],tmp[nn])) #tmp[nn] last element is tag_noRec = tot-Riak...
 } }}}  #loop iak ends 
 
 Rs <- as.data.frame(Rs); names(Rs) <- c('i','a','k','g','j','k1','R','phi','delta')  # this is simulated rec data
 
 print(dim(Rs))
 
 
 all(Rs$R>0)
 Ts <- as.data.frame(Ts); names(Ts) <- c('i','a','k','tag_tot','tag_noRec')  # this is simulated tag data (basically the observed)
 
 sum(Ts[,ncol(Ts)]>0)
 sum(E$tag_noRec>0)
     
 # now match aggregated rec to tag
 # match rec+cat data to tag data by iak
 iak_tag <- paste(Ts$i,Ts$a,Ts$k,sep='_')
 length(iak_tag)==length(unique(iak_tag))  #T
 all(iak_tag==E$tag_iak) #T
 iak_rec <- paste(Rs$i,Rs$a,Rs$k,sep='_') 
 tmp <- aggregate(Rs$R, list(iak_rec), 'sum'); names(tmp) <- c('iak','R') #aggregating rec over iak groups
  all(tmp$iak == levels(as.factor(iak_rec)))
  # convert levels of rec IAK to integers
  rec_iak <- as.integer(as.factor(iak_rec))
  all(tmp$iak %in% iak_tag) #T
  !all(iak_tag %in% tmp$iak)  #T
  findLoc <- function(target, refs){
   inds <- which(refs==target)
   if(length(inds)==0) return(NA)
   if(length(inds)==1) return(inds)
   if(length(inds)>1) stop('Multuple match')
  }
  ind_rec2tag <- as.integer(sapply(iak_tag, findLoc, tmp$iak))
  head(cbind(Ts, tmp[ind_rec2tag,]),10)
  Ts$R <- tmp$R[ind_rec2tag];  Ts$R[is.na(Ts$R)] <- 0
  ind_rec2tag[is.na(ind_rec2tag)] <- nrow(tmp) + 1 #a trick: replace missing index (tagged iak, but not found in rec) with n (number of iak found) + 1. Next, add 0 to the value. So value[index] will do 2 things: match back, replace missing values with 0.
  all(Ts$R == c(tmp$R,0)[ind_rec2tag]) # see how this trick works
  
  tag_noRec <- Ts$tag_noRec
  save(file='true',Rs=Rs,Ts=Ts,Ps=Ps,ind_rec_KJAG=ind_rec_KJAG,tag_noRec=tag_noRec,ind_rec2tag=ind_rec2tag)
  
  
  
  # vectorization to avoid for loops
  rec_lags <- Rs$j - Rs$i 
 # vectorize multiple indices for rec+cat data to avoid for loop
ind_rec_KJAG <- ind_rec_KJA <- ind_rec_JK <- rep(NA, nrow(Rs)) #ind_rec_KKAJ <- : I believe we no longer need it
for(i0 in 1:nrow(Rs)){
  i <- Rs$i[i0]
  j <- Rs$j[i0]
  k <- Rs$k[i0]
  k1 <- Rs$k1[i0]
  a <- Rs$a[i0]
  a1 <- min( a+ rec_lags[i0], A )  #Rs$agep[i0]
  g <- Rs$g[i0]
  
  i1 <- 1 #if(j<=4) i1 <- 1 else i1 <- 2
  if(k1<=3) i2 <- 1 else if(k1<=5) i2 <- 2 else if(k1<=8) i2 <- 3 else if(k1==9) i2 <- 4 else i2 <- 5
  tmp1 <- matrix(0,J2,K2); tmp1[i1,i2] <- 1
  ind_rec_JK[i0] <- which(tmp1==1)  
  tmp1 <- array(0, c(K1,J,A,G)); tmp1[k1,j,a1,g] <- 1
  ind_rec_KJAG[i0] <- which(tmp1==1)
  tmp1 <- array(0, c(K1,J,A)); tmp1[k1,j,a1] <- 1
  ind_rec_KJA[i0] <- which(tmp1==1) 
#  tmp1 <- array(0, c(K,K1,A,J)); tmp1[k,k1,a1,j] <- 1
#  ind_rec_KKAJ[i0] <- which(tmp1==1)
}

  #********************  overwrite session 2 ******************
  E$tag_tot <- Ts$tag_tot
  E$tag_noRec <- Ts$tag_noRec
  E$tag_iak <- iak_tag
  E$ind_rec2tag <- ind_rec2tag
  E$nrec <- nrow(Rs)
  E$rec_R <- Rs$R
  E$rec_IJKK1AG <- Rs[,c('i','j','k','k1','a','g')]
  E$rec_phi_ij <- Rs$phi
  E$rec_delta_gjk <- Rs$delta
  E$rec_iak <- rec_iak
  E$rec_lags <- Rs$j - Rs$i
  E$ind_rec_KJAG <- ind_rec_KJAG
  E$ind_rec_KJA <- ind_rec_KJA
  E$ind_rec_JK <- ind_rec_JK
  #************************************************************
  
  # calculate R_tot, R_m at rec data level. note this is not actually observed, we only observed sum of them for each gjk, not iakgjk.
  Rs$lambda <- lambda_jk[ind_rec_JK]
  Rs$R_m <- rep(NA,nrow(Rs)); for(i in 1:nrow(Rs)) Rs$R_m[i] <- rbinom(1,Rs$R[i],Rs$delta[i]/(Rs$delta[i]+(1-Rs$delta[i])*Rs$lambda[i]))
  #Rs$delta/(Rs$delta+(1-Rs$delta)*Rs$lambda)
  Rs$R_um <- Rs$R-Rs$R_m 
                                        
  gjk_cat <- paste(E$cat_KJG$gear,E$cat_KJG$year,E$cat_KJG$basin,sep='_')
  length(gjk_cat)==length(unique(gjk_cat))
  gjk_rec <- paste(Rs$g,Rs$j,Rs$k,sep='_')
  cat_R_tot <- cat_R_um <- rep(NA,length(catch_m)) # recall this is the actually observed catch_m
  for(i in 1:length(catch_m)){
    inds <- which(gjk_rec==gjk_cat[i])
    cat_R_tot[i] <- sum(Rs$R[inds])
    cat_R_um[i] <- sum(Rs$R_um[inds])
  }   
  Cs$R_tot_simu <- cat_R_tot
  Cs$R_um_simu <- cat_R_um
  
  lambda <- lambda_jk[E$ind_cat_JK]
  delta <- E$cat_delta_gjk[E$ind_cat_positive]      
  
  # vectorization for cat data
  tmp <- Cs[Cs$R_tot_simu>0 & Cs$catch_m>0, ]
  ind_cat_JK <- rep(NA, nrow(tmp))
  for(i in 1:nrow(tmp)){
    i1 <- 1 #if(tmp$year[i]<=4) i1 <- 1 else i1 <- 2
    if(tmp$k1[i]<=3) i2 <- 1 else if(tmp$k1[i]<=5) i2 <- 2 else   if(tmp$k1[i]<=8) i2 <- 3 else if(tmp$k1[i]==9) i2 <- 4 else i2 <- 5
    tmp1 <- matrix(0,J2,K2); tmp1[i1,i2] <- 1
    ind_cat_JK[i] <- which(tmp1==1)  
  }
  
  #********************  overwrite session 3 ******************
  E$cat_R_tot <- cat_R_tot
  E$cat_R_um <- cat_R_um
  E$ind_cat_JK <- ind_cat_JK
  E$ind_cat_positive <- which(cat_R_tot>0 & catch_m>0)
  #************************************************************  
  
  #print(loglik) 
  
return(E)
}# function ends 
  
  
  
  