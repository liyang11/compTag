

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
E$lbs <- x$lbs; E$ubs <- x$ubs; E$flags <- as.integer(x$flags); E$nams <- x$nams; E$dims <- x$dims;
x <- x$paras

#x1 <- apply(b,2,mean)
#b <- t(rbind(apply(b,2,mean),apply(b,2,quantile,c(0.025,0.975))))
#mat1 <- cbind(as.numeric(a$x0), b)
#
#write.csv(file='out.csv',mat1,row.names=F)
#
#source('translate.R')
##x <- as.numeric(a$x0)
#b <- translate(x1, E, model=1)
#b$M_jak
#
##b[E$flags==5,]


source('translate.R')
out <- as.list(rep(NA,6)); k0 <- 1
for(m in 1:2){
a <- readMat(paste('out',m,'_0.mat',sep=''))
b <-  a$matpara
#plot(b[,which(E$flags==5)], type='l')
b <- t(rbind(apply(b,2,mean),apply(b,2,quantile,c(0.025,0.975))))
nams <- paste('Model',m, nams0 <- c('mean','lower','upper'),sep='_')
for(i in 1:length(nams0)){
   tmp <- translate(b[,i], E, model=m)
   out[[k0]] <- tmp; names(out)[k0] <- nams[i]; k0 <- k0+1
 }
}




# print posterior mean estimates
foo <- function(saveTxt=FALSE){

for(m in 1:2){
if(saveTxt) sink(paste('out',m,'.txt',sep=''))
cat('\n\n************ Posterior mean (lower and upper 95% credible interval) of parameter estimates for model',m,' ***********\n')
for(j in 1:9){
cat('\n',names(out[[1]])[j],':\n',sep='')
for(i in 1:length(nams0)){
cat(nams0[i],':\n',sep='')
tmp <- as.data.frame(out[[(m-1)*3+i]][[j]]); if(length(tmp)==1) names(tmp) <- NULL; printRowNam <- TRUE; if(nrow(tmp)==1) printRowNam <- FALSE
print(tmp, row.names=printRowNam); cat('\n')
}
cat('\n')
}
print(out[[(m-1)*3+1]][10:length(out[[1]])])
cat('************************************************\n\n\n')
if(saveTxt) sink()
}

}

foo(saveTxt=TRUE)






#m <- readMat('paras1.mat')$mat
#p <- dim(m)[1]
#nsim <- dim(m)[3]-1
#Cs <- matrix(NA,p,nsim); for(i in 1:nsim) Cs[,i] <- as.integer(m[,3,i+1]<=m[,1,i+1]&m[,4,i+1]>=m[,1,i+1])
##tmp <- order(apply(Cs,2,mean),decreasing=T)
#cp <- apply(Cs,1,mean)
#Ms <- apply(m[,2,-1],1,mean)
#
#opt_q <- 2 #catchability: 1=q, 2=qj, 3=qk, 4=qjk
#opt_M <- 1 #natural mortality: 1=M, 2=Mj, 3=Ma, 4=Mk 
#opt_Pi <- 1 #dynamic: 1=Pi, 2=Pi_a
#parametrizations <- c(opt_q,opt_M,opt_Pi)
#names(parametrizations) <- c('q(1=q,2=qj,3=qk)','M(1=M,2=Mj,3=Ma,4=Mk)','Pi(1=Pi,2=Pi_a)')
#
#source('readData.R')  # read data into environmental var. E
#E <- readData()
#E$parametrizations <- parametrizations
#
#source('initPara.R')  # initialize parameter 
#x <- initPara(E)
##lapply(x,head,20)
##write.csv(file='x.csv',x[c(5,1:4)],row.names=F)
#
## we purify x to be the parameter vector and copy its other information to environmental variable E
#E$lbs <- x$lbs; E$ubs <- x$ubs; E$flags <- x$flags; E$nams <- x$nams; E$dims <- x$dims;
#
#a <- as.data.frame(cbind(m[,1,2],Ms,cp)); row.names(a) <- E$nams 
#write.csv(file='a.csv',a)
#
#y <- m[,2,1]
#y[E$flags==5]
#
#m0 <- m[,,1]
#names()





