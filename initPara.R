
initPara <- function(E){
#require('gtools') #for R function: rdirichlet
#model_id <- as.numeric(commandArgs(TRUE))
#model_id <- 1

opt_q <- E$parametrizations[1] #1=q, 2=qj, 3=qk, 4=qjk
opt_M <- E$parametrizations[2] #1=M, 2=Mj, 3=Ma, 4=Mk 
opt_Pi <- E$parametrizations[3] #1=Pi, 2=Pi_a

J <- E$J; K <- E$K; A <- E$A; K1 <- E$K1; K0 <- E$K0

#============================ initialize parameters
set.seed(20)

nams <- nams0 <- character(); flags <- integer(); paras <- lbs <- ubs <- numeric();

# initial Lambda using empirical estimates from cat data
inds <- which(E$cat_R_tot>0)
tmp <- E$cat_KJG[inds,]
tmp$R_m <- E$cat_R_tot[inds] - E$cat_R_um[inds]
tmp$R_um <- E$cat_R_um[inds]
tmp$catch <- exp(E$log_cat_catch[inds])
tmp$catch_m <- tmp$catch*E$cat_delta_gjk[inds]
uid <- as.factor(paste(tmp$year, tmp$basin,sep='_'))
lid <- length(levels(uid))
#x0 <- tmp$R_um*tmp$delta_gjk/(tmp$R_m*(1-tmp$delta_gjk))
R_m <- R_um <- array(0, c(E$J2, E$K2, 2))
for(i in 1:nrow(tmp)){
  i1 <- 1 #if(tmp$year[i]<=4) i1 <- 1 else i1 <- 2
  if(tmp$basin[i]<=3) i2 <- 1 else if(tmp$basin[i]<=5) i2 <- 2 else if(tmp$basin[i]<=8) i2 <- 3 else if(tmp$basin[i]==9) i2 <- 4 else i2 <- 5 
  R_m[i1,i2,1] <- R_m[i1,i2,1] + tmp$R_m[i]
  R_m[i1,i2,2] <- R_m[i1,i2,2] + tmp$catch_m[i]
  R_um[i1,i2,1] <- R_um[i1,i2,1] + tmp$R_um[i]
  R_um[i1,i2,2] <- R_um[i1,i2,2] + tmp$catch[i]  
}
Delta <- R_m[,,2]/R_um[,,2]
lambda_jk <- R_um[,,1]*Delta/(R_m[,,1]*(1-Delta))
lambda_jk[lambda_jk>1] <- 0.9
lambda_jk[lambda_jk==0] <- 0.1
nams0 <- c(nams0, 'lambda')
nams <- c(nams, paste('lambda[', 1:length(lambda_jk),']',sep=''))
flags <- c(flags, rep(1,length(lambda_jk))); paras <- c(paras, lambda_jk); lbs <- c(lbs, rep(0,length(lambda_jk))); ubs <- c(ubs, rep(1,length(lambda_jk)))


# load initial guess
#source('int_N_q.R')
load('Inits1')
a1 <- as.data.frame(E$rec_levels$rec_basin); names(a1) <- 'basin'; a1$y <- 1:nrow(a1)
a1 <- merge(a1, C_k)
N_1k <- numeric(nrow(a1)); N_1k[a1$y] <- a1$N_1k
N_1ak_mu <- rowMeans(N_1ak)
N_1ak_dev <- N_1ak; for(i in 1:ncol(N_1ak)) N_1ak_dev[,i] = N_1ak_dev[,i]/N_1ak_mu

if(opt_q==1) {q_jk <- -9; nams <- c(nams,'q')}
if(opt_q==2) {q_jk <- rep(-9,J); nams <- c(nams, paste('q[j=',1:J,']',sep=''))}
if(opt_q==3) {q_jk <- rep(-9,K); q_jk <- rep(-9,K); nams <- c(nams, paste('q[k=',1:K1,']',sep=''))}
# 2016-8-25: reparametrize q_jk instead of a scalar
#q_jk <- log(q_jgk) #log(0.0000015) #-13.41  #log q ~ Uniform(-15, -10)
#q_jk <- rep(log(q_jgk), J*K); 
if(opt_q==4) {q_jk <- rep(-9, J*K); nams <- c(nams, paste(rep('q',J*K),'[j=',rep(1:J,K),',k=',rep(1:K,each=J),']',sep=''))}
nams0 <- c(nams0, 'q')

flags <- c(flags, rep(2,length(q_jk)))
paras <- c(paras, q_jk)
lbs <- c(lbs, rep(-15,length(q_jk)))
ubs <- c(ubs, rep(-8,length(q_jk)))

tau <- 1.26 #tau ~ Uniform(-3, 3) 
eta <- log(13.074) #this is log eta ~ Uniform(-3, 5)
nams0 <- c(nams0, 'tau','eta')
nams <- c(nams, 'tau', 'eta') 
flags <- c(flags, 3, 4) 
paras <- c(paras, tau, eta) 
lbs <- c(lbs, c(-1,-1)) 
ubs <- c(ubs, c(4,5))


#M_jak <- log(0.25)  #M_jak <- array(0.25, dims[c(6, 5, 2, 4)]) #KJAG

# 2016-8-25: reparametrize M_jak as a function of age
agevec <- E$agevec
nams0 <- c(nams0, 'M')
#M_jak <- log(0.25)/mean(agevec)  # this is log(x) so that M=x^a
## 2016-9-8: I would try M_jak = M_a
m0 <- 0.25
if(opt_M==1) M_jak <- log(m0)
if(opt_M==2) M_jak <- rep(log(m0),J)
if(opt_M==3) M_jak <- rep(log(m0),A) #M_jak*agevec
#M_jak <- matrix(log(0.25), J, length(agevec)) 
if(opt_M==4) M_jak <- rep(log(m0),K1) #M_jak*agevec
tmp <- c('j','a','k')
if(opt_M==1) nams <- c(nams,'M')
if(opt_M>1) nams <- c(nams, paste('M[',tmp[opt_M-1],'=',1:length(M_jak),']',sep=''))
flags <- c(flags, rep(5,length(M_jak)))
paras <- c(paras, as.vector(M_jak))
lbs <- c(lbs,rep(-5,length(M_jak)))
ubs <- c(ubs, rep(2,length(M_jak))) 


# we work on unconstrained H instead of Pi
if(opt_Pi==1) Pi0 <- array(0, c(K,K1,1))
if(opt_Pi==2) Pi0 <- array(0, c(K,K1,A))
H0 <- array(0, c(K,K1-1,dim(Pi0)[3]))
for(a in 1:dim(Pi0)[3]){
 for(k in 1:K) {Pi0[k,k,a] <- .9; Pi0[k,-k,a] <- .1/(K-1)}
 H0[,,a] <- matrix(0, K, K1-1)
 for(k in 1:K){H0[k,,a] <- log(Pi0[k,1:(K1-1),a]/Pi0[k,K1,a])}
}


nams0 <- c(nams0, 'H(Pi)')
if(opt_Pi==1) nams <- c(nams, paste('H[k=', rep(1:K, K1-1), ',k1=',rep(1:(K1-1),each=K),']', sep=''))
if(opt_Pi==2) nams <- c(nams, paste('H[k=', rep(rep(1:K, K1-1),A), ',k1=',rep(rep(1:(K1-1),each=K),A),',a=',rep(1:A,each=K*(K1-1)),']', sep=''))
flags <- c(flags, rep(6,length(H0)))
paras <- c(paras, as.vector(H0))
lbs <- c(lbs, rep(-30,length(H0))) #Inf
ubs <- c(ubs, rep(30,length(H0))) #Inf


#mu1 <- log(N_1k) #rep((17.5+8)/2, K) #runif(K, 8.1, 13) #log(N), a=1
#dev1 <- matrix(0, K, J)  #a=1 runif(K*J, -7.8, 7.9)
#mu2 <- log(N_1ak_mu) #rep((16+8)/2, K) #runif(K, 8.1, 13) #log(N), j=1
#dev2 <- log(N_1ak_dev) #matrix(0, K, A-1) #j=1  runif(K*(A-1), -7.8, 7.9)
nu0 <- 0.1


# 2016-6-15: CAUTION: use 1:K0 only
tmp <- matrix(rep(N_1k, J), K, J)
R_kj <- log(tmp[1:K0,])  #recruitment for all tagging years 
N0_ak <- log(N_1ak[1:K0,]) #initial abundance at age in the first tagging year

nams0 <- c(nams0, 'R','N0','nu0')
nams <- c(nams, 
 paste('R[k=', rep(1:K0, J), ',j=',rep(1:J,each=K0),']', sep=''),
 paste('N0[k=', rep(1:K0, A-1), ',a=',rep(2:A,each=K0),']', sep=''),'nu')
flags <- c(flags, rep(7,length(R_kj)), rep(8,length(N0_ak)), rep(9, length(nu0)) )
paras <- c(paras, as.vector(R_kj), as.vector(N0_ak), as.vector(nu0))
lbs <- c(lbs, rep(-20,length(R_kj)), rep(-20,length(N0_ak)), rep(exp(-5), length(nu0)) ) 
ubs <- c(ubs, rep(log(6e11),length(R_kj)), rep(log(5e9),length(N0_ak)), rep(exp(5), length(nu0)) )   
                                
dims <- table(flags)
names(dims) <- nams0

return(list(paras=paras,lbs=lbs,ubs=ubs,flags=flags,nams=nams,dims=dims))

} # end of function
