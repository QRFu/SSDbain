library(SSDbain)
library(testthat)
library(stringr)
library(conflicted)
library(pkgmaker, warn.conflicts = FALSE)
library(foreach)
library(doParallel)
library(doRNG)
library(pkgmaker, warn.conflicts = FALSE)
no_cores <- detectCores() - 1
## Loading required package: iterators
cl <- makeCluster(no_cores)

registerDoParallel(cl)

clusterEvalQ(cl, .libPaths())



#========================================================================================
#  using SSDbain to test the BF02 for the two-sided Welch's t-test if H0 is true
#========================================================================================

Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=1000)
N_J1<-Res[[1]][[5]]
eta1_J1<-Res[[1]][[3]]
eta2_J1<-Res[[1]][[4]]

N_J2<-Res[[2]][[5]]
eta1_J2<-Res[[2]][[3]]
eta2_J2<-Res[[2]][[4]]


N_J3<-Res[[3]][[5]]
eta1_J3<-Res[[3]][[3]]
eta2_J3<-Res[[3]][[4]]
#sample size required when J<-1
#Using N=87 and b

#========================================================================================
#  J=1
#========================================================================================
J<-1
#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H0 is true
#========================================================================================")
m1<-0
m2<-0
var1<-1
var2<-1
sd1<-sqrt(var1)
sd2<-sqrt(var2)
BFthresh<-3
N_SIM<-1000
Variances<-'equal'
N1<-N2<-N_J1
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)


bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF02_J1<-k/N_SIM



#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H2 is true
#========================================================================================")
m1<-0.5
m2<-0

#BF20<-calBF02medium1(m1,m2,sd1,sd2,N,J,BFthresh,N_SIM,Variances,eta)
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)

bf<-1/bf

bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF20_J1<-k/N_SIM




#========================================================================================
#  J=2
#========================================================================================
J<-2
#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H0 is true
#========================================================================================")
m1<-0
m2<-0
BFthresh<-3
Variances<-'equal'
N1<-N2<-N_J2
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)


bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF02_J2<-k/N_SIM



#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H2 is true
#========================================================================================")
m1<-0.5
m2<-0

#BF20<-calBF02medium1(m1,m2,sd1,sd2,N,J,BFthresh,N_SIM,Variances,eta)
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)

bf<-1/bf

bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF20_J2<-k/N_SIM


#========================================================================================
#  J=3
#========================================================================================
J<-3
#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H0 is true
#========================================================================================")
m1<-0
m2<-0
BFthresh<-3
Variances<-'equal'
N1<-N2<-N_J3
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)


bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF02_J3<-k/N_SIM




#========================================================================================")
#  using SSDbain to test the BF20 for the two-sided Welch's t-test if H2 is true
#========================================================================================")
m1<-0.5
m2<-0

#BF20<-calBF02medium1(m1,m2,sd1,sd2,N,J,BFthresh,N_SIM,Variances,eta)
samph<-c(N1,N2)
ngroup<-samph/J

##simulate data
##Repeat the same simulation results
set.seed(10)
##simulate data
datalist<-list()
for(i in 1:N_SIM){
  datalist[[i]]<-cbind(rnorm(N1,m1,sd1),rnorm(N2,m2,sd2))
}

##caculate estimate and variance for group mean
estimate<-lapply(datalist, function(x) apply(x,2,mean))
variance<-lapply(datalist, function(x) apply(x,2,function(x) as.matrix(var(x)/N1) ) )
##caculate bf
bf<-c()

BF<- foreach(i = 1:N_SIM) %dorng% {
  library(bain)

  #  for(i in 1:N_SIM){
  names(estimate[[i]]) <- c("mu1","mu2")
  variance_sp<-(variance[[i]][1]+variance[[i]][2])/2
  covlist<-list(matrix(variance_sp),matrix(variance_sp))
  res<-bain(estimate[[i]], n=ngroup, "mu1 = mu2",Sigma=covlist,group_parameters=1,joint_parameters = 0)
  bf<-res$fit$BF[1]
  return(bf)
}
bf<-unlist(BF)

bf<-1/bf

bf<-bf[!is.na(bf)]
N_SIM<-length(bf)
k<-0
for(i in 1:N_SIM){
  if (bf[[i]]>BFthresh)
    k<-k+1
}
P_BF20_J3<-k/N_SIM


test_that("SSDbain", {expect_equal(eta1_J1, P_BF02_J1,tolerance = .01)})
test_that("SSDbain", {expect_equal(eta2_J1, P_BF20_J1,tolerance = .01)})


test_that("SSDbain", {expect_equal(eta1_J2, P_BF02_J2,tolerance = .01)})
test_that("SSDbain", {expect_equal(eta2_J2, P_BF20_J2,tolerance = .01)})


test_that("SSDbain", {expect_equal(eta1_J3, P_BF02_J3,tolerance = .01)})
test_that("SSDbain", {expect_equal(eta2_J3, P_BF20_J3,tolerance = .01)})


