
rm(list=ls())

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
tmp <- E; tmp <- tmp[-which(names(tmp)=='rec_levels')]
writeMat(paste('Drq',opt_q,'m',opt_M,'p',opt_Pi,'.mat',sep=''), E=tmp, x=x, L0=L0)




#*********  save simulation data
nsim <- 1
for(repl in 1:nsim){ #******** simulate data starts
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

} #******** simulate data ends

cputime <- as.numeric(proc.time()[3]-ptm)
cputime <- cputime/60
cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')






simpleSummary <- FALSE
if(simpleSummary){
  repl <- 1
  #x1 <- readMat('out1.mat')$xtmp
  x2 <- readMat('paras1.mat')
  #mat <- as.data.frame(cbind(x,as.vector(x1),x2,E$lbs,E$ubs))
  mat <- as.data.frame(cbind(x2$mat,E$lbs,E$ubs))
  row.names(mat) <- E$nams
  names(mat) <- c('initial','estB','estB_l','estB_u','lbs','ubs') #'estF',
  write.csv(file='paraF.csv',mat)
  print('done')
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









