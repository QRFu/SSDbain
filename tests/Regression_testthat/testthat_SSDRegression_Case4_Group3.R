rm(list=ls())
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
#========================================================================================")
#  using SSDanova to test the BF01 for the two-sided ANOVA if H0 is true
#========================================================================================")
N_min<-10
N_max<-1000
T_sim<-1000
#R_square<-0.13#effect size
k<-3#number of predictors
l=0.2
rho<-matrix(runif(k*k),nrow=k)
for (i in 1:k){
  for (j in 1:k){
    if(i==j){
      rho[i,j]=1
    }else
      rho[i,j]=l
  }
}
#rho<-matrix(c(1,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2,0.2,1),nrow=4)
Hyp1<-'beta1>0&beta2>0&beta3>0'
Hyp2<-'Hc'
ratio<-c(1,1,1)
R_square1=0.13
R_square2=0.13
BFthresh=3
eta=0.8
seed=10
standardize=FALSE
Res<-SSDRegression(Hyp1,Hyp2,k,rho,R_square1,R_square2,T_sim,BFthresh,eta,seed,standardize,ratio)
N_J1<-Res[[1]][1]
eta1_J1<-Res[[2]][[1]][[1]]
eta2_J1<-Res[[2]][[1]][[2]]


T<-T_sim
#========================================================================================
#  J=1
#========================================================================================
J<-1
N<-N_J1
beta<-cal_beta(k,R_square1,R_square2,Hyp1,Hyp2,ratio,rho)

beta1<-beta[[1]]
beta2<-beta[[2]]
varnames<-c()
b<-character()
for (i in 1:k){
  varnames[i]=paste('beta',as.character(i),sep='')
}

num_eq <- sum(unlist(gregexpr("[=]",Hyp1))>0)
num_0 <- sum(unlist(gregexpr("[0]",Hyp1))>0)
num_larger <- sum(unlist(gregexpr("[>]",Hyp1))>0)
#bf0j_H0<-cal_bf(N,k,beta1,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
beta<-beta1
set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion")
#beta<-c(d,d,d)
bf<-numeric(T)
samph<-rep(N,k)
sampm<-samph/J
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
flag_Hyp = 1

BF<- foreach(ii = 1:T, .packages = c("bain")) %dorng% {
  #  for (ii in 1:T){
  res<- list()
  library(MASS)
  X<-mvrnorm(n=N,mu=numeric(k),Sigma=rho)
  if(flag_Hyp){
    sd<-sqrt(1-R_square1)
  }else{
    sd<-sqrt(1-R_square2)
  }
  epsilon<-rnorm(n=N,mean=0,sd=sd)
  #beta1<-beta2<-beta3<-d
  Y<-  epsilon
  for (i in 1:k)
  {
    Y<-Y+beta[i]*X[,i]
  }
  library(fungible)
  library(stats)


  if(standardize==TRUE){
    # Obtain standardized regression coefficients
    capture.output(int <-fungible:: seBeta(X = X, y = Y, Nobs = N,
                                           estimator = "Normal"), file='NUL')
    # Store the parameter estimates
    estimate <- int$CIs[, 2]
    for (i in 1:k)
    {
      names(estimate)[i]<-paste('beta',i,sep='')
    }
    #names(estimate)<-c('beta1','beta2','beta3')
    # Store the variance covariance matrix for the parameters
    # of interest
    cov <- list(int$cov.mat)
    cov<-matrix(unlist(cov),nrow = k)
  }else{
    M=cbind(X,Y)
    data<-as.list(as.data.frame(M))
    for (i in 1:k)
    {
      names(data)[i]<-paste('beta',i,sep='')
    }
    Equ<-'Y ~ beta1'
    for (i in 2:k){
      Equ<-paste(Equ,' + beta',i,sep='')
    }
    Equ<-paste('lm(',Equ,'-1',',data)',sep='')
    # execute a multiple regression using lm()
    regr<-eval(parse(text=Equ))
    # take a look at the estimated regression coefficients and their names
    coef(regr)
    estimate <- coef(regr)
    cov<-vcov(regr)
  }

  if(flag_Hyp){

    capture.output(res<-bain::bain(estimate,Hyp1,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')
  }else{
    capture.output(res<-bain::bain(estimate,Hyp2,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')

  }

  if(Hyp2=='Hc'){
    bf<-res$fit$BF[[1]]
  } else{
    bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  return(bf)
}
bf<-unlist(BF)
bf0u_H0<-bf

bf0j_H0<-na.omit(bf0u_H0)
T0<-length(bf0j_H0)
medbf0j_H0<-median(bf0j_H0)
kk<-lapply(1:T0, function (x) bf0j_H0[x]>BFthresh )
eta1_test_J1<-sum(unlist(kk)*1)/T0


flag_Hyp = 1
beta<-beta2
set.seed(seed,kind="Mersenne-Twister",normal.kind="Inversion")
#beta<-c(d,d,d)
bf<-numeric(T)
samph<-rep(N,k)
sampm<-samph/J
BF<- foreach(ii = 1:T, .packages = c("bain")) %dorng% {
  #  for (ii in 1:T){
  res<- list()
  library(MASS)
  X<-mvrnorm(n=N,mu=numeric(k),Sigma=rho)
  if(flag_Hyp){
    sd<-sqrt(1-R_square1)
  }else{
    sd<-sqrt(1-R_square2)
  }
  epsilon<-rnorm(n=N,mean=0,sd=sd)
  #beta1<-beta2<-beta3<-d
  Y<-  epsilon
  for (i in 1:k)
  {
    Y<-Y+beta[i]*X[,i]
  }
  library(fungible)
  library(stats)


  if(standardize==TRUE){
    # Obtain standardized regression coefficients
    capture.output(int <-fungible:: seBeta(X = X, y = Y, Nobs = N,
                                           estimator = "Normal"), file='NUL')
    # Store the parameter estimates
    estimate <- int$CIs[, 2]
    for (i in 1:k)
    {
      names(estimate)[i]<-paste('beta',i,sep='')
    }
    #names(estimate)<-c('beta1','beta2','beta3')
    # Store the variance covariance matrix for the parameters
    # of interest
    cov <- list(int$cov.mat)
    cov<-matrix(unlist(cov),nrow = k)
  }else{
    M=cbind(X,Y)
    data<-as.list(as.data.frame(M))
    for (i in 1:k)
    {
      names(data)[i]<-paste('beta',i,sep='')
    }
    Equ<-'Y ~ beta1'
    for (i in 2:k){
      Equ<-paste(Equ,' + beta',i,sep='')
    }
    Equ<-paste('lm(',Equ,'-1',',data)',sep='')
    # execute a multiple regression using lm()
    regr<-eval(parse(text=Equ))
    # take a look at the estimated regression coefficients and their names
    coef(regr)
    estimate <- coef(regr)
    cov<-vcov(regr)
  }

  if(flag_Hyp){

    capture.output(res<-bain::bain(estimate,Hyp1,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')
  }else{
    capture.output(res<-bain::bain(estimate,Hyp2,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')

  }

  if(Hyp2=='Hc'){
    bf<-res$fit$BF[[1]]
  } else{
    bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
  }
  return(bf)
}
bf0u_Hj<-unlist(BF)


bfj0_Hj<-1/bf0u_Hj

bfj0_Hj<-na.omit(bfj0_Hj)
T0<-length(bfj0_Hj)
# medianbfj0_Hj<-median(bfj0_Hj)
kk<-lapply(1:T0, function (x) bfj0_Hj[x]>BFthresh )
eta2_test_J1<-sum(unlist(kk))/T0



# #test Tables in paper
test_that("SSDbain", {expect_equal(eta1_test_J1,eta1_J1 ,tolerance = .01)})
test_that("SSDbain", {expect_equal(eta2_test_J1,eta2_J1,tolerance = .01)})

