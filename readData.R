
readData <- function(){  
  
# ********************** Part 1
# data import and organization 
dat <- as.list(rep(NA,4))
dat[[1]] <- read.csv('tag_iak.csv')
dat[[2]] <- read.csv('rec_iakgjk.csv')
dat[[3]] <- read.csv('cat_eff_rep_new.csv') #'cat_eff_rep.csv'
dat[[4]] <- read.csv('age_gjka.csv')

names(dat) <- c('tag','rec','cat','age')
nd <- length(dat)
ng <- c(3,6,3,4) #last column of grouping variables
minAge <- min(as.integer(dat$age$agep))

# adjustment. CAUTION!
dat$rec[c(236, 244),]  #small and large gill
dat$rec$R[236] <- 2
dat$rec <- dat$rec[-244, ]

# data processing # 2016-6-21: exclude two rec_basins
excludeSites <- TRUE
if(excludeSites){
 table(dat[[2]]$rec_basin)
 dat[[2]] <- dat[[2]][-which(dat[[2]]$rec_basin %in% c('Georgia_bay','North_channel')),]
 dat[[2]]$rec_basin <- droplevels(dat[[2]]$rec_basin)
 dat[[3]] <- dat[[3]][-which(dat[[3]]$basin %in% c('Georgia_bay','North_channel')),]
 dat[[3]]$basin <- droplevels(dat[[3]]$basin)
 dat[[4]] <- dat[[4]][-which(dat[[4]]$rec_basin %in% c('Georgia_bay','North_channel')),]
 dat[[4]]$rec_basin <- droplevels(dat[[4]]$rec_basin)
}

lapply(dat, dim)
nrec <- nrow(dat$rec)
ncat <- nrow(dat$cat) 

## data processing # 2016-9-20: exclude rec_Y=2011
#dat[[2]] <- dat[[2]][-which(dat[[2]]$rec_Y==2011), ]
#dat[[3]] <- dat[[3]][-which(dat[[3]]$year==2011), ]  
#dat[[4]] <- dat[[4]][-which(dat[[4]]$year==2011), ]     


# check data
all(dat$rec$rec_Y >= dat$rec$tag_Y)

# align levels
# CAUTION: the levels of tag year and recovery year must have overlapping level (same starting year). Need check
lev1 <- levels(dat$tag$tag_basin)
lev2 <- levels(dat$cat$basin)
lev3 <- c(lev2[lev2 %in% lev1], lev2[!(lev2 %in% lev1)])

dat$tag$tag_basin <- factor(as.character(dat$tag$tag_basin), levels=lev3)
dat$rec$tag_basin <- factor(as.character(dat$rec$tag_basin), levels=lev3)
dat$rec$rec_basin <- factor(as.character(dat$rec$rec_basin), levels=lev3)
dat$cat$basin <- factor(as.character(dat$cat$basin), levels=lev3)
dat$age$rec_basin <- factor(as.character(dat$age$rec_basin), levels=lev3)

levels(dat$rec$rec_gear) <- c('Gill','Gill','Trap')
levels(dat$age$gear) <- c('Gill','Trap')

# convert types, integerize the group variables
lev <- as.list(rep(NA,nd)); names(lev) <- names(dat)
for(i in 1:nd){
 lev[[i]] <- as.list(rep(NA,ng[i])) 
 for(j in 1:ng[i]){
  if(!is.factor(dat[[i]][,j])) 
    dat[[i]][,j] <- as.factor(dat[[i]][,j])
  lev[[i]][[j]] <- levels(dat[[i]][,j])
  names(lev[[i]]) <- names(dat[[i]])[1:ng[i]]
  dat[[i]][,j] <- as.integer(dat[[i]][,j])
 }
}

#all(dat$rec$rec_Y >= dat$rec$tag_Y)

lapply(dat, 'head')

dat$tag$iak <- paste(dat$tag[,1],dat$tag[,2],dat$tag[,3],sep='_')
#all(table(dat$tag$iak)==1)
dat$rec$iak <- paste(dat$rec[,1],dat$rec[,2],dat$rec[,3],sep='_')
dat$rec$gjk <- paste(dat$rec[,4],dat$rec[,5],dat$rec[,6],sep='_')
#all(table(paste(dat$rec$iak, dat$rec$gjk, sep='-'))==1)

all(table(paste(dat$rec$iak, dat$rec$gjk, sep='-'))==1)

# precalculation
agevec <- 4:12
phiFun <- function(x) return(1-0.31/(1+exp(8.52-x)))
rec_phi_ij <- phiFun( (dat$rec$rec_Y - dat$rec$tag_Y)*12 + 1 )

dat$cat$delta_gjk <- dat$cat$catch_m/dat$cat$catch
#all(table(paste(dat$cat$gear, dat$cat$year, dat$cat$basin, sep='_'))==1)  

dims <- as.integer(unlist(lapply(lev$rec, length)))
#I A K G J K
#1 2 3 4 5 6 
I <- dims[1]
A <- dims[2]
K <- dims[3]
G <- dims[4]
J <- dims[5]
K1 <- dims[6]
J2 <- 1 #2 #groups of J
K2 <- 5; if(excludeSites) K2 <- 3 #groups of K1
K0 <- 8 



findLoc <- function(target, refs){
 inds <- which(refs==target)
 if(length(inds)==0) return(NA)
 if(length(inds)==1) return(inds)
 if(length(inds)>1) stop('Multuple match')
}
  
# match cat data to rec data by gjk          
cat_gjk <- paste(dat$cat$gear, dat$cat$year, dat$cat$basin, sep='_')
all(dat$rec$gjk%in%cat_gjk)#T
!all(cat_gjk%in%dat$rec$gjk)#T
ind_cat2rec <- sapply(dat$rec$gjk, findLoc, cat_gjk)
all(dat$rec$gjk==cat_gjk[ind_cat2rec])  #T
#-------- get delta_gjk for rec data from delta_gjk in cat data
rec_delta_gjk <- dat$cat$delta_gjk[ind_cat2rec] #for dat$rec

## a useful check: 
#tmp <- rep(NA, nrec); for(i in 1:nrec) tmp[i] <- sum(dat$rec$R[dat$rec$gjk==dat$rec$gjk[i]])
##head(dat$cat[ind_cat2rec,])
##head(dat$rec)
##sum(dat$rec$R[dat$rec$gjk=='2_1_2'])
##head(tmp)
#all(dat$cat$R_tot[ind_cat2rec]==tmp)

# CAUTION: after merge, foo2 reordered dat$rec!!!
#foo2 <- merge(x=dat$rec, y=foo2, by='gjk', all=F)
rec_lags <- dat$rec$rec_Y - dat$rec$tag_Y

# match rec+cat data to tag data by iak  
tmp <- aggregate(dat$rec$R, list(dat$rec$iak), 'sum'); names(tmp) <- c('iak','R') #aggregating rec over iak groups
all(tmp$iak == levels(as.factor(dat$rec$iak)))
# convert levels of rec IAK to integers
rec_iak <- as.integer(as.factor(dat$rec$iak))
all(tmp$iak %in% dat$tag$iak) #T
!all(dat$tag$iak %in% tmp$iak)  #T
tag_rec <- dat$tag # to add some variables such as total catch for each iak, I don't want to change dat$tag
length(tag_rec$iak)==length(unique(tag_rec$iak))  #T
ind_rec2tag <- as.integer(sapply(tag_rec$iak, findLoc, tmp$iak))
head(ind_rec2tag,10)  #if tag iak not found in aggregated rec iak (tmp), then it will be NA
head(cbind(dat$tag, tmp[ind_rec2tag,]),10)
head(tag_rec$iak)
tag_rec$R <- tmp$R[ind_rec2tag]
tag_rec$R[is.na(tag_rec$R)] <- 0 
sum(tmp$R[tmp$iak=='1_3_1'])
tag_noRec <- tag_rec$T-tag_rec$R #this is T_iak - R_iak...
tag_iak <- tag_rec$iak # this is only needed for simulation
#foo <- merge(x=tmp, y=dat$tag, by='iak', all=T)
#foo$R[is.na(foo$R)] <- 0
#diffs <- foo$T-foo$R

#ind_match <- rep(NA, nrow(tag_rec)); for(i in 1:nrow(tag_rec)) if(any(tmp$iak == tag_rec$iak[i])) ind_match[i] <- which(tmp$iak == tag_rec$iak[i])
ind_rec2tag[is.na(ind_rec2tag)] <- nrow(tmp) + 1 #a trick: replace missing index (tagged iak, but not found in rec) with n (number of iak found) + 1. Next, add 0 to the value. So value[index] will do 2 things: match back, replace missing values with 0.
all(tag_rec$R == c(tmp$R,0)[ind_rec2tag]) # see how this trick works

cat_KJG <- dat$cat[,c('basin','year','gear')]
log_cat_catch <- log(dat$cat$catch)
cat_delta_gjk <- dat$cat$delta_gjk
#sum(abs(cat_catch*cat_delta_gjk-dat$cat$catch_m))
cat_effort <- dat$cat$effort
cat_R_tot <- dat$cat$R_tot 
cat_R_um <- dat$cat$R_um
#all(cat_R_tot-cat_R_um==dat$cat$R_m)

#logEpart <- dat$age$n_sum*dat$age$p_perc
# 2016-7-14
#logEpart <- 25*dat$age$p_perc 
age_n_gjka <- 25*dat$age$p_perc #old logEpart




# **************************** Part 2
# vectorize multiple indices for rec+cat data to avoid for loop
ind_rec_KJAG <- ind_rec_KJA <- ind_rec_JK <- rep(NA, nrow(dat$rec)) #ind_rec_KKAJ <- : I believe we no longer need it
for(i0 in 1:nrow(dat$rec)){
  i <- dat$rec$tag_Y[i0]
  j <- dat$rec$rec_Y[i0]
  k <- dat$rec$tag_basin[i0]
  k1 <- dat$rec$rec_basin[i0]
  a <- dat$rec$agep[i0]
  a1 <- min( dat$rec$agep[i0]+ rec_lags[i0], A )
  g <- dat$rec$rec_gear[i0]
  
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
#ind_rec_KKAJ <- ind_rec_KKAJ[which(lags==0)]
rec_IJKK1AG <- dat$rec[,c('tag_Y','rec_Y','tag_basin','rec_basin','agep','rec_gear')]


# vectorize j,k for cat data to directly access lambda_jk
#tmp <- dat$cat[dat$cat$R_tot>0, ]
tmp <- dat$cat[dat$cat$R_tot>0 & dat$cat$catch_m>0, ]
ind_cat_JK <- rep(NA, nrow(tmp))
for(i in 1:nrow(tmp)){
  i1 <- 1 #if(tmp$year[i]<=4) i1 <- 1 else i1 <- 2
  if(tmp$basin[i]<=3) i2 <- 1 else if(tmp$basin[i]<=5) i2 <- 2 else   if(tmp$basin[i]<=8) i2 <- 3 else if(tmp$basin[i]==9) i2 <- 4 else i2 <- 5
  tmp1 <- matrix(0,J2,K2); tmp1[i1,i2] <- 1
  ind_cat_JK[i] <- which(tmp1==1)  
}
ind_cat_positive <- which(dat$cat$R_tot>0 & dat$cat$catch_m>0)

# Avoid for loop to grealy speed up
ind_age_KJGA <- ind_age_KJG <- rep(NA, nrow(dat$age))
for(i0 in 1:nrow(dat$age)){
 k1 <- dat$age$rec_basin[i0]
 j <- dat$age$year[i0]
 g <- dat$age$gear[i0]
 a <- dat$age$agep[i0]
 tmp <- array(0, c(K,J,G,A)); tmp[k1,j,g,a] <- 1
 ind_age_KJGA[i0] <- which(tmp==1)
 tmp <- array(0, c(K,J,G)); tmp[k1,j,g] <- 1
 ind_age_KJG[i0] <- which(tmp==1)
}

#induse_age <- which(dat$age$rec_basin <= K0)  # no need since now K0 == K1

# for debug
tag_dbIAK <- tag_rec[,1:3] #tag_Y,agep,tag_basin 

# save environmental variable
E <- list(
 I=I, A=A, K=K, J=J, G=G, K1=K1, J2=J2, K2=K2, K0=K0,
 nrec=nrec, ncat=ncat, agevec=as.numeric(agevec),
 tag_iak=tag_iak, tag_noRec=tag_noRec, ind_rec2tag=ind_rec2tag,
 tag_dbIAK=tag_dbIAK, 
 tag_tot=tag_rec$T, #only need for simulation 
 rec_levels=lev$rec, #only need for initials using load('Inits1')
 rec_R=as.numeric(dat$rec$R), rec_IJKK1AG=rec_IJKK1AG, rec_phi_ij=rec_phi_ij, rec_delta_gjk=as.numeric(rec_delta_gjk), rec_iak=rec_iak, rec_lags=as.numeric(rec_lags), ind_rec_KJAG=ind_rec_KJAG, ind_rec_KJA=ind_rec_KJA, ind_rec_JK=ind_rec_JK, 
 cat_KJG=cat_KJG, log_cat_catch=log_cat_catch, cat_delta_gjk=as.numeric(cat_delta_gjk), cat_effort=as.numeric(cat_effort), cat_R_tot=as.numeric(cat_R_tot), cat_R_um=as.numeric(cat_R_um), ind_cat_positive=ind_cat_positive, ind_cat_JK=ind_cat_JK, 
 age_n_gjka=as.numeric(age_n_gjka), ind_age_KJGA=ind_age_KJGA, ind_age_KJG=ind_age_KJG
)

return(E)
} # end of function

## save for checking
#tmp <- foo; tmp$diffs <- diffs
#tmp$ind_match <- ind_match
#write.csv(file='tag_comb.csv',tmp, row.names=F)
#table(dat$rec$tag_basin,dat$rec$rec_basin)
#tmp <- foo2; tmp$ind_rec_2 <- ind_rec_2
#write.csv(file='rec_comb.csv',tmp, row.names=F)
#tmp <- dat$cat[dat$cat$R_tot>0 & dat$cat$catch_m>0, ]
#tmp$ind_cat_2 <- ind_cat_2
#tmp <- tmp[order(ind_cat_2),]
#write.csv(file='cat_comb.csv',tmp, row.names=F)





