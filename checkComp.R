
rm(list=ls())

# choose parametrizations for the model to fit
opt_q <- 2 #catchability: 1=q, 2=qj, 3=qk, 4=qjk
opt_M <- 1 #natural mortality: 1=M, 2=Mj, 3=Ma, 4=Mk 
opt_Pi <- 1 #dynamic: 1=Pi, 2=Pi_a
parametrizations <- c(opt_q,opt_M,opt_Pi)
names(parametrizations) <- c('q(1=q,2=qj,3=qk)','M(1=M,2=Mj,3=Ma,4=Mk)','Pi(1=Pi,2=Pi_a)')

source('readData.R')  # read data into environmental var. E
E <- readData()
E$parametrizations <- parametrizations

source('initPara.R')  # initialize parameter 
x <- initPara(E)

# we purify x to be the parameter vector and copy its other information to environmental variable E
E$lbs <- x$lbs; E$ubs <- x$ubs; E$flags <- x$flags; E$nams <- x$nams; E$dims <- x$dims;
x <- x$paras



repl <- 2
source('simuData.R')
seeds <- repl*10
E1 <- simuData(x,E,seeds)
              

source('getLogLik.R')
a1 <- getLogLik(x,E1,comp=1:4,model=1)
as.character(a1)
-sum(a1)

source('getLogLik.R')
a2 <- getLogLik(x,E1,comp=1:4,model=2)
as.character(a2)
-sum(a2)

#load('F1'); F1=Fs
#load('F2'); F2=Fs
#all(F1==F2)

#rm(list=ls())
a <- as.list(rep(NA,3))
load('true')
a[[1]] <- list(Ps=Ps,Rs=Rs,Ts=Ts,ind_rec_KJAG=ind_rec_KJAG,tag_noRec=tag_noRec,ind_rec2tag=ind_rec2tag)

tmp <- aggregate(a[[1]]$Ps, list(E1$rec_iak), 'sum'); names(tmp) <- c('iak','P')
tmpP <- c(tmp$P,0)[E1$ind_rec2tag]

sum(a[[1]]$Rs[,7]*log(a[[1]]$Ps)) + sum(E1$tag_noRec*log(1-tmpP))

load('m1')
a[[2]] <- list(Ps=Ps,Rs=Rs,Ts=Ts)
tmp <- aggregate(a[[2]]$Ps, list(E1$rec_iak), 'sum'); names(tmp) <- c('iak','P')
tmpP <- c(tmp$P,0)[E1$ind_rec2tag]
sum(a[[2]]$Rs[,1]*log(a[[2]]$Ps)) + sum(E1$tag_noRec*log(1-tmpP))

load('m2')
a[[3]] <- list(Ps=Ps,Rs=Rs,Ts=Ts,ind_rec_KJAG=ind_rec_KJAG)
tmp <- aggregate(a[[3]]$Ps, list(E1$rec_iak), 'sum'); names(tmp) <- c('iak','P')
tmpP <- c(tmp$P,0)[E1$ind_rec2tag]
sum(a[[3]]$Rs[,1]*log(a[[3]]$Ps)) + sum(E1$tag_noRec*log(1-tmpP))


a1 <- a[[1]]$Rs[,7]*log(a[[1]]$Ps)
a2 <- a[[2]]$Rs[,1]*log(a[[2]]$Ps)
a3 <- a[[3]]$Rs[,1]*log(a[[3]]$Ps)
sum(a1)
sum(a2)
sum(a3)

all(a1 == a3)

mat <- cbind(a1, a2)  #, a[[2]]$Ps
write.csv(file='foo.csv',mat,row.names=F)
