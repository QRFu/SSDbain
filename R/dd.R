cal_bf<-function(N,k,beta,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp,seed,standardize){
  #   N1<-N2<-N3<-N
  #simulate data
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
    # M=cbind(X,Y)
    # data<-as.list(as.data.frame(M))
    # for (i in 1:k)
    # {
    #   names(data)[i]<-paste('beta',i,sep='')}
    # # load the bain package which includes the simulated sesamesim data set
    # library(bain)
    #
    # Equ<-'Y ~ beta1'
    # for (i in 2:k){
    #   Equ<-paste(Equ,' + beta',i,sep='')
    # }
    # Equ<-paste('lm(',Equ,',data)',sep='')
    # # execute a multiple regression using lm()
    # regr<-eval(parse(text=Equ))
    # test hypotheses with bain. Note that standardized = FALSE denotes that the
    # hypotheses are in terms of unstandardized regression coefficients
    if(flag_Hyp){
    #res<-bain::bain(regr,n=N/J,'beta1 = 0 & beta2=0 & beta3=0',standardize = TRUE)
   # res<-bain::bain(regr,'beta1 = 0 & beta2=0 & beta3=0',fraction = 1,standardize = TRUE)
    #res<-bain::bain(regr,Hyp1,fraction = 1,standardize = FALSE)
       #capture.output(res<-bain::bain(regr, hypothesis=Hyp1,fraction = J, standardize = standardize), file='NUL')
     # capture.output(res<-bain::bain(regr,hypothesis= Hyp1,fraction = J, standardize = standardize), file='NUL')

       capture.output(res<-bain::bain(estimate,Hyp1,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')

    }else{
      #capture.output(res<-bain::bain(regr,hypothesis= Hyp2,fraction = J, standardize = standardize), file='NUL')
      capture.output(res<-bain::bain(estimate,Hyp2,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')

    }

    #capture.output(res<-bain::bain(estimate,Hyp2,fraction = J,n=N,Sigma=cov,group_parameters=0,joint_parameters=k), file='NUL')

     # res<-bain::bain(regr,'beta1>0 & beta2>0 & beta3>0',fraction = 2, standardize = TRUE)

    #  res<-bain::bain(regr,n=N/J,'beta1>0 & beta2>0 & beta3>0',standardize = TRUE)
    #res<-bain::bain(estimate,Hyp2,n=sampm,Sigma=Sigma,group_parameters=0,joint_parameters = k)

    #bf1<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    if(Hyp2=='Hc'){
      bf<-res$fit$BF[[1]]
    } else{
      bf<-res$fit$Fit[[1]]/res$fit$Com[[1]]
    }
     return(bf)
  }
  bf<-unlist(BF)
   return(bf)
}

