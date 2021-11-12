
SSDRegression<-function(Hyp1,Hyp2,k,rho,R_square1,R_square2,T_sim=1000,BFthresh,eta,seed=10,standardize=TRUE,ratio){
  cat('It may take 30-45 minutes to compute the required sample size. Please wait...')

  if(BFthresh<1){
    stop("BFthresh should be larger than 1!")
  }
  if(eta<0){
    stop("eta should be larger than 0")
  }
  sym<-unlist(strsplit(Hyp1,split='[<]'))
  len<-length(sym)
  if(len>1){
    stop("Please check the Hyp1!",call. = TRUE)
  }
  if(Hyp2=='Ha'){
     if(R_square1!=0){
       stop("R_square1 does not match with the Hyp1!",call. = TRUE)
     }
  }
if(sum(abs(ratio))>sum(ratio)){
  stop("Please check the hypothesis and the ratio!",call. = TRUE)
}

  if(length(R_square1)==k&&length(R_square2)==k){
    beta1<-R_square1
    beta2<-R_square2
    beta<-numeric(k)
    # beta=seq(from=k,to=1,by=-1)
    beta_r1<-matrix(runif(k*k),nrow = k)
    beta_r2<-matrix(runif(k*k),nrow = k)
    for (i in 1:k){
      for (j in 1:k){
        beta_r1[i,j]<-beta1[i]*beta1[j]*rho[i,j]
        beta_r2[i,j]<-beta2[i]*beta2[j]*rho[i,j]
      }
    }
    R_square1=sum(sum(beta_r1))
    R_square2=sum(sum(beta_r2))

  }else if(length(R_square1)==1&&length(R_square2)==1){
  beta<-cal_beta(k,R_square1,R_square2,Hyp1,Hyp2,ratio,rho)
  beta1<-beta[[1]]
  beta2<-beta[[2]]
  }else{
  cat('Please input correct value for R_square or beta')
  }
  N_min<-10
  N_max<-1000
  N_J<-c()
  medBF_J<-list()
  N_min_J<-c(N_min,N_min,N_min)
  N_max_J<-c(N_max,N_max,N_max)
  options(warn=-1)
  # library(parallel)
  library(stringr)
  library(conflicted)
  library(pkgmaker, warn.conflicts = FALSE)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(pkgmaker, warn.conflicts = FALSE)
  no_cores <- detectCores() - 2
  ## Loading required package: iterators
  cl <- makeCluster(no_cores)

  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths())
  varnames<-c()
  b<-character()
  for (i in 1:k){
    varnames[i]=paste('beta',as.character(i),sep='')
  }

  num_eq <- sum(unlist(gregexpr("[=]",Hyp1))>0)
  num_0 <- sum(unlist(gregexpr("[0]",Hyp1))>0)
  num_larger <- sum(unlist(gregexpr("[>]",Hyp1))>0)
  if(num_eq>=k-1){
    Jmax<-3}else{
    Jmax<-1
  }
  for (J in 1:Jmax){
    if(J>1){
      cat(paste('###Sensitive Analysis  ',as.character(J),'b###','\n',sep=''))
    }
    else
    {      cat(paste('###Sensitive Analysis  ', 'b###','\n',sep=''))
    }
    ##step1 sample size for the population that the sample exactly equal to  the population mean and variance
    #increase sample size step by step till the BF reach the threshold
    N_min<-N_min_J[J]
    N_max<-N_max_J[J]
    # cat(paste('Initial Step: fast determine sample size with the real mean and variance.','\n',sep=''))
    # capture.output(Results_Step1<-Fast_SSD(N_min,N_max,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b), file='NUL')
    # N_Step1<-Results_Step1[[1]]
    # medBF0j_H0_Step1<-Results_Step1[[2]]
    # medBFj0_Hj_Step1<-Results_Step1[[3]]
    # cat(paste('The initial sample size is ',as.character(N_Step1),'\n',sep=''))
    # cat(paste('Median BF0j is ',as.character(signif(medBF0j_H0_Step1,4)),'\n',sep=''))
    # cat(paste('Median BFj0 is ',as.character(signif(medBFj0_Hj_Step1,4)),'\n','\n',sep=''))
    #
    # #step2 substitute the N_Step1 to the function
    # cat(paste('Second Step: exactly determine sample size.','\n',sep=''))
    N_Step1<-100
    medBF<- cal_medbf(N_Step1,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

    #medBF<- cal_medbf(N_Step1,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
    medBF0j_H0_Step2<-medBF[[1]]
    medBFj0_Hj_Step2<-medBF[[2]]
    cat(paste('The sample size is ',as.character(N_Step1),'\n',sep=''))
    cat(paste('P(BF0j>',as.character(BFthresh),') is ',as.character(signif(medBF[[1]],4)),'\n',sep=''))
    cat(paste('P(BFj0>',as.character(BFthresh),') is ',as.character(signif(medBF[[2]],4)),'\n','\n',sep=''))

    if(medBF[[1]]>eta&&medBF[[2]]>eta)
    {
      N_max<-N_Step1
      #cat(paste('The median BF is larger than the threshold. Therefore,','\n',sep=''))
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
      N_min<-max(floor(N_max/2),10);
      medBF<- cal_medbf(N_min,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

     # medBF<- cal_medbf(N_min,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      N_temp<-c()
      while((medBF[[1]]>eta&&medBF[[2]]>eta)&&N_min>10){
        N_temp<-N_min
        N_min<-max(floor(N_min/2),10);
        medBF<- cal_medbf(N_min,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)
        #medBF<- cal_medbf(N_min,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
        cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n',sep=''))
      }
      if(N_min*2<(N_Step1-1))
      {
        N_max<-N_temp
      }
      else
        N_max<-N_Step1#min(N_min*2,N_Step1)
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))
    }
    else{
      N_max<-N_Step1
      cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n','\n',sep=''))

      while(!(medBF[[1]]>eta&&medBF[[2]]>eta)){
        N_max<-2*N_max;
        medBF<- cal_medbf(N_max,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

        #medBF<- cal_medbf(N_max,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
        cat(paste('The upper boundary of sample size is N_max=',as.character(N_max),'\n',sep=''))
      }
      if(N_max>N_Step1*2+1)
        N_min<-N_max/2
      else
        N_min<-N_Step1#max(N_max/2,N_Step1)
      cat(paste('The lower boundary of sample size is N_min=',as.character(N_min),'\n','\n',sep=''))
    }
    flag_end<-0
    if(N_min==10)
    {
      medBF<- cal_medbf(N_min,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

      #medBF<- cal_medbf(N_min,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
      if(medBF[[1]]>eta&&medBF[[2]]>eta)
      {
        flag_end<-1
        N<-N_min
        medBF_N10<-medBF
        #cat(paste('N=',as.character(N),'\n',sep=''))
      }
    }


    N_med<-floor((N_min+N_max)/2)
    while((N_med-N_min>=1&&N_max-N_med>=1)&&!flag_end){
      cat(paste('The sample size in current loop is N=',as.character(N_med),'\n',sep=''))
      medBF<- cal_medbf(N_med,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

      #medBF<- cal_medbf(N_med,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
      cat(paste('P(BF0j>',as.character(BFthresh),') is ',as.character(signif(medBF[[1]],4)),'\n',sep=''))
      cat(paste('P(BFj0>',as.character(BFthresh),') is ',as.character(signif(medBF[[2]],4)),'\n','\n',sep=''))

      if(medBF[[1]]>eta&&medBF[[2]]>eta){
        N_max<-N_med
      }
      else{
        N_min<-N_med
      }
      N_med<-floor((N_min+N_max)/2)
    }
    if(flag_end){
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<-medBF_N10

    }
    else{
      N<-N_max
      cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))
      medBF<- cal_medbf(N,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

      #medBF<- cal_medbf(N,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)

      # medBF_Vec<-c(medBF[[1]],medBF[[2]])
      # E_Vec<-c(medBF[[3]],medBF[[4]])
      # WE<-(medBF[[9]]+medBF[[10]])/2
      # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
      # CI_Vec<-c(medBF[[7]],medBF[[8]])
      # names(medBF_Vec)<-c('medBF01','medBF10')
      # names(E_Vec)<-c('Error I','Error II')
      # names(M_Vec)<-c('ME I ','ME II','WE')
      # print(medBF_Vec)
      # cat(paste('\n',sep=''))
      # print(E_Vec)
      # cat(paste('\n',sep=''))
      # print(M_Vec)
      # cat(paste('\n',sep=''))
      # cat(paste('CI1',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[7]])
      # cat(paste('CI2',sep=''))
      # cat(paste('\n',sep=''))
      # print(medBF[[8]])

    }

    # medBF<- cal_medbf(N_med,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
    # if(medBF[[1]]>Med_b&&medBF[[2]]>Med_b){
    #   N<-N_med
    # }
    # else{
    #   N<-N_max
    #   cat(paste('N=',as.character(N),'\n',sep=''))
    #   medBF<- cal_medbf(N,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
    #
    # }
    N_J[J]<-N
    medBF_J[[J]]<-medBF
  }
  for (J in 1:Jmax){
    N<-N_J[J]
    medBF<-medBF_J[[J]]
    medBF_Vec<-c(medBF[[1]],medBF[[2]])
    # E_Vec<-c(medBF[[3]],medBF[[4]])
    # WE<-(medBF[[9]]+medBF[[10]])/2
    # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
    # CI_Vec<-c(medBF[[7]],medBF[[8]])
    if(Jmax>1 && J>1)
      cat(paste('#####Sensitive analysis ',as.character(J), '*', 'b','#####','\n',sep=''))

    cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))

    names(medBF_Vec)<-c('eta01','eta10')
    # names(E_Vec)<-c('Error I','Error II')
    # names(M_Vec)<-c('ME I ','ME II','WE')
    print(medBF_Vec)
    cat(paste('\n',sep=''))
    # print(E_Vec)
    # cat(paste('\n',sep=''))
    # print(M_Vec)
    # cat(paste('\n',sep=''))
    # cat(paste('quantile1',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[7]])
    # cat(paste('quantile2',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[8]])

  }
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('Hypothesis:','\n',sep=''))
  cat(paste(Hyp1,' vs ',Hyp2,'\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('BFthreth:','\n',sep=''))
  cat(paste(as.character(BFthresh),'\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('eta:','\n',sep=''))
  cat(paste(as.character(eta),'\n',sep=''))

  stopCluster(cl)

  return(list(N_J,medBF_J))

}



  CalRegression<-function(N,Hyp1,Hyp2,k,rho,R_square1,R_square2,T_sim,BFthresh,eta,seed=10,standardize=TRUE){
  # Hyp1<-'beta1 = 0 & beta2=0 & beta3=0'
  # Hyp2<-'beta1 > 0 & beta2>0 & beta3>0'
    if(length(R_square1)==k&&length(R_square2)==k){
      beta1<-R_square1
      beta2<-R_square2
    }else if(length(R_square1)==1&&length(R_square2)==1){
      beta<-cal_beta(k,R_square1,R_square2,Hyp1,Hyp2)
      beta1<-beta[[1]]
      beta2<-beta[[2]]
    }else{
      cat('Please input correct value for R_square or beta')
    }
  N_J<-c()
  medBF_J<-list()
  N_min_J<-c(N_min,N_min,N_min)
  N_max_J<-c(N_max,N_max,N_max)
  options(warn=-1)
  N_min<-10
  N_max<-1000
  # library(parallel)
  library(stringr)
  library(conflicted)
  library(pkgmaker, warn.conflicts = FALSE)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(pkgmaker, warn.conflicts = FALSE)
  no_cores <- detectCores() - 2
  ## Loading required package: iterators
  cl <- makeCluster(no_cores)

  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths())
  Jmax<-3
  for (J in 1:Jmax){
    if(J>1){
      cat(paste('###Sensitive Analysis  ',as.character(J),'b###','\n',sep=''))
    }
    else
    {      cat(paste('###Sensitive Analysis  ', 'b###','\n',sep=''))
}
    ##step1 sample size for the population that the sample exactly equal to  the population mean and variance
    #increase sample size step by step till the BF reach the threshold
    N_min<-N_min_J[J]
    N_max<-N_max_J[J]
    # cat(paste('Initial Step: fast determine sample size with the real mean and variance.','\n',sep=''))
    # capture.output(Results_Step1<-Fast_SSD(N_min,N_max,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b), file='NUL')
    # N_Step1<-Results_Step1[[1]]
    # medBF0j_H0_Step1<-Results_Step1[[2]]
    # medBFj0_Hj_Step1<-Results_Step1[[3]]
    # cat(paste('The initial sample size is ',as.character(N_Step1),'\n',sep=''))
    # cat(paste('Median BF0j is ',as.character(signif(medBF0j_H0_Step1,4)),'\n',sep=''))
    # cat(paste('Median BFj0 is ',as.character(signif(medBFj0_Hj_Step1,4)),'\n','\n',sep=''))
    #
    # #step2 substitute the N_Step1 to the function
    # cat(paste('Second Step: exactly determine sample size.','\n',sep=''))
    N_Step1<-N
    medBF<- cal_medbf(N_Step1,k,beta1,beta2,R_square1,R_square2,T_sim,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize)

    #medBF<- cal_medbf(N_Step1,mu1,mu2,sd,T_sim,J,varnames,hyp1,hyp2,ERr2,IRr2,Med_b,flag_fast,type)
    medBF0j_H0_Step2<-medBF[[1]]
    medBFj0_Hj_Step2<-medBF[[2]]
    cat(paste('The sample size is ',as.character(N_Step1),'\n',sep=''))
    cat(paste('P(BF0j>',as.character(BFthresh),') is ',as.character(signif(medBF[[1]],4)),'\n',sep=''))
    cat(paste('P(BFj0>',as.character(BFthresh),') is ',as.character(signif(medBF[[2]],4)),'\n','\n',sep=''))

    # }
    N_J[J]<-N
    medBF_J[[J]]<-medBF
  }
  for (J in 1:Jmax){
    N<-N_J[J]
    medBF<-medBF_J[[J]]
    medBF_Vec<-c(medBF[[1]],medBF[[2]])
    # E_Vec<-c(medBF[[3]],medBF[[4]])
    # WE<-(medBF[[9]]+medBF[[10]])/2
    # M_Vec<-c(medBF[[5]],medBF[[6]],WE)
    # CI_Vec<-c(medBF[[7]],medBF[[8]])
    if(Jmax>1 && J>1)
      cat(paste('#####Sensitive analysis ',as.character(J), '*', 'b','#####','\n',sep=''))

    cat(paste('The final sample size is N=',as.character(N),'\n',sep=''))

    names(medBF_Vec)<-c('eta01','eta10')
    # names(E_Vec)<-c('Error I','Error II')
    # names(M_Vec)<-c('ME I ','ME II','WE')
    print(medBF_Vec)
    cat(paste('\n',sep=''))
    # print(E_Vec)
    # cat(paste('\n',sep=''))
    # print(M_Vec)
    # cat(paste('\n',sep=''))
    # cat(paste('quantile1',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[7]])
    # cat(paste('quantile2',sep=''))
    # cat(paste('\n',sep=''))
    # print(medBF[[8]])

  }
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('Hypothesis:','\n',sep=''))
  cat(paste(Hyp1,' vs ',Hyp2,'\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('BFthresh:','\n',sep=''))
  cat(paste(as.character(BFthresh),'\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('\n',sep=''))
  cat(paste('eta:','\n',sep=''))
  cat(paste(as.character(eta),'\n',sep=''))



  stopCluster(cl)

  return(list(N_J,medBF_J))

}
