
#rm(list=ls())

repl <- as.numeric(commandArgs(TRUE))
print(repl)

require(R.matlab)
require(truncnorm) 

#repl <- as.numeric(commandArgs(TRUE))
ptm <- proc.time()[3]
  
# choose parametrizations for the model to fit
opt_q <- 1 #catchability: 1=q, 2=qj, 3=qk, 4=qjk
opt_M <- 1 #natural mortality: 1=M, 2=Mj, 3=Ma, 4=Mk 
opt_Pi <- 1 #dynamic: 1=Pi, 2=Pi_a
parametrizations <- c(opt_q,opt_M,opt_Pi)
names(parametrizations) <- c('q(1=q,2=qj,3=qk)','M(1=M,2=Mj,3=Ma,4=Mk)','Pi(1=Pi,2=Pi_a)')

source('readData.R')  # read data into environmental var. E
E <- readData()
E$parametrizations <- parametrizations

lapply(E,head, 20)





source('initPara.R')  # initialize parameter 
x <- initPara(E)
#lapply(x,head,20)
#write.csv(file='x.csv',x[c(5,1:4)],row.names=F)

# we purify x to be the parameter vector and copy its other information to environmental variable E
E$lbs <- x$lbs; E$ubs <- x$ubs; E$flags <- x$flags; E$nams <- x$nams; E$dims <- x$dims;
x <- x$paras

x[E$flags==5]

#lapply(E,class)

source('getLogLik.R')
a <- getLogLik(x,E,comp=1:4,model=1)
as.character(a)
L0 <- sum(a)
##-sum(a)

#*********  save real data
#tmp <- E; tmp <- tmp[-which(names(tmp)=='rec_levels')]
#writeMat(paste('Drq',opt_q,'m',opt_M,'p',opt_Pi,'.mat',sep=''), E=tmp, x=x, L0=L0)




#*********  save simulation data
nsim <- 50
#for(repl in 1:nsim){ #******** simulate data starts
cat(repl,' '); if(!repl%%20) cat('\n')
 
# simulation
source('simuData.R')
seeds <- repl*10
E1 <- simuData(x,E,seeds)


source('getLogLik.R')
a1 <- getLogLik(x,E1,comp=1:4,model=1)
as.character(a1)
sum(a1)

source('getLogLik.R')
a2 <- getLogLik(x,E1,comp=1:4,model=2)
as.character(a2)
sum(a2)

L1 <- sum(a1)

tmp <- E1; tmp <- tmp[-which(names(tmp)=='rec_levels')]
writeMat(paste('Ds',repl,'q',opt_q,'m',opt_M,'p',opt_Pi,'.mat',sep=''), E=tmp, x=x, L0=L1)

#} #******** simulate data ends

cputime <- as.numeric(proc.time()[3]-ptm)
cputime <- cputime/60
cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')






simpleSummary <- FALSE  # summary for simulation
if(simpleSummary){
 require(R.matlab)
 nsim <- 50
 for(m in 1:2){
  for(i in 1:nsim){
   a <- readMat(paste('out',m,'_',i,'.mat',sep=''))
   b <-  a$matpara; b <- t(rbind(apply(b,2,mean),apply(b,2,quantile,c(0.025,0.975))))
   if(i==1&&m==1) tmp <- t(a$x0)
   if(i==1) ms <- rs <- cp <- matrix(0,nsim,np <- length(a$x0))
   ms[i,] <- b[,1]
   rs[i,] <- (b[,1]-a$x0)^2
   cp[i,] <- b[,2]<a$x0 & b[,3]>a$x0
  }
  sds <- apply(ms,2,'sd')
  ms <- colMeans(ms) # posterior mean averaging over nsim replications
  rs <- sqrt(colMeans(rs)) # RMSE
  cp <- colMeans(cp) # Coverage proability
  tmp <- cbind(tmp,ms,sds,rs,cp)
 }
 
 # get nams
 opt_q <- 1 #catchability: 1=q, 2=qj, 3=qk, 4=qjk
 opt_M <- 1 #natural mortality: 1=M, 2=Mj, 3=Ma, 4=Mk 
 opt_Pi <- 1 #dynamic: 1=Pi, 2=Pi_a
 parametrizations <- c(opt_q,opt_M,opt_Pi)
 names(parametrizations) <- c('q(1=q,2=qj,3=qk)','M(1=M,2=Mj,3=Ma,4=Mk)','Pi(1=Pi,2=Pi_a)')
 source('readData.R')  # read data into environmental var. E
 E <- readData()
 E$parametrizations <- parametrizations
 source('initPara.R')  # initialize parameter 
 x <- initPara(E)
 nams <- x$nams 
  
 tmp <- as.data.frame(tmp); names(tmp) <- c("True", paste("Model",rep(c(1,2),each=4), '_',
  rep(c('Mean','std','RMSE','CP'),2), sep='')); row.names(tmp) <- nams
 write.csv(file='outS.csv',tmp)
 mean(tmp$Model2_CP >= tmp$Model1_CP)
 mean(tmp$Model2_RMSE <= tmp$Model1_RMSE)
}







#getLogLik(x,E1,comp=1:4)
#fit <- optim(paras, getLogLik, model=1, method='L-BFGS-B',lower=lbs,upper=ubs,control=list(trace=1, REPORT=1,factr=1e7,maxit=5),hessian=F)


#x <- readMat('tmp.mat')$x2
##as.character(lik <- getLogLik(x, comp=1:3))
##-sum(lik)
##
##x[flags==7]
##paras[flags==7]
##
##getLogLik(x, comp=1:3)
##getLogLik(paras, comp=1:3)
##
###getLogLik(paras, comp=c(2))
###t0 <- proc.time()
##
##
##
###for(i in 1:4)  getLogLik(paras, comp=c(3,4))
###proc.time()-t0
##
##t0 <- proc.time()
##fit <- optim(paras, getLogLik, model=model_id, method='L-BFGS-B',lower=lbs,upper=ubs,control=list(trace=1, REPORT=1,factr=1e7,maxit=5),hessian=F) # tolerance of about 1e-8.
##proc.time()-t0
##
##source('translate.R')
##x <- translate(fit$par)
##
###translate(paras)
##
##save(file=paste('out',model_id,sep=''),x=x, fit=fit)
##
##
###cbind(paras,lbs,ubs)
##              
##              
###inds <- which(flags<6)
###paras2 <- paras
###paras2[inds] <- log((paras2[inds]-lbs[inds])/(ubs[inds]-paras2[inds]))
###source('getLogLik2.R')
###
###fit <- optim(paras2, getLogLik2, method='Nelder-Mead', control=list(trace=1, REPORT=1,reltol=1e-16, maxit=10),hessian=T)
##
##
###sum(abs(paras[inds] - (lbs[inds] + (ubs[inds]-lbs[inds])/(1+exp(-paras2[inds])) )))
##
##
###logit <- function(x,a,b) return(log((x-a)/(b-x)))
###vec <- seq(1e-4, 2-1e-4,len=1000)
###plot(vec, logit(vec,1e-4, 2-1e-4), type='b')
##
##
###len <- length(x0 <- as.vector(lambda_jk))
###getL2(x0)
###
####x0 <- rep(0.5, len)
####getL2(x0)
###
###fit <- optim(x0, getL2, method='L-BFGS-B',lower=rep(0.001,len),upper=rep(0.999,len),control=list(trace=1)) #,control=list(trace=0,reltol=1e-16, maxit=1e3)
###
###fit$par
##
##
##
### read results from MATLAB
##obj <- readMat('out.mat')
##
##source('translate.R')
##
##x <- obj$x
##exp(x[flags==5])
##
##
##-getLogLik(obj$x, comp=1:3)
##
##x[flags==5] <- 0.2
##-getLogLik(x, comp=1:3)
##
##
##writeMat('tmp.mat',x=x)
##
##y <- translate(obj$x)









