cal_medbf<-function(N,k,beta1,beta2,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,BFthresh,seed,standardize){
  if(Hyp2=='Hc'){
    bf0j_H0<-cal_bf(N,k,beta1,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
  }else if(Hyp2=='Ha'){
    bf0j_H0<-cal_bf(N,k,beta1,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
  }else{
    #BF0u(BFiu)
    bf0u_H0<-cal_bf(N,k,beta1,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
    bfju_H0<-cal_bf(N,k,beta1,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 0,seed,standardize)
    bf0j_H0<-bf0u_H0/bfju_H0
  }
  bf0j_H0<-na.omit(bf0j_H0)
  T0<-length(bf0j_H0)
  medbf0j_H0<-median(bf0j_H0)
  kk<-lapply(1:T0, function (x) bf0j_H0[x]>BFthresh )
   eta1<-sum(unlist(kk)*1)/T0
  # kk<-lapply(1:T0, function (x) bf0j_H0[x]<1/3 )
  # M_H0<-sum(unlist(kk))/T0
  # kk<-lapply(1:T0, function (x) bf0j_H0[x]>1/3&bf0j_H0[x]<3 )
  # W_H0<-sum(unlist(kk))/T0
  # CI_H0<-round(quantile(bf0j_H0,c(0.05,0.1,0.2,0.8,0.9,0.95)),2)

  #cat(paste('Median BF0j = ',as.character(medbf0j_H0),'\n',sep=''))


  #Hj is true
  if(Hyp2=='Hc'){
    bf0j_Hj<-cal_bf(N,k,beta2,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
    bfj0_Hj<-1/bf0j_Hj
  }else if(Hyp2=='Ha'){
    bf0j_Hj<-cal_bf(N,k,beta2,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
    bfj0_Hj<-1/bf0j_Hj
    } else{
  #BF0u(BFiu)
    bf0u_Hj<-cal_bf(N,k,beta2,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 1,seed,standardize)
    bfju_Hj<-cal_bf(N,k,beta2,R_square1,R_square2,T,rho,J,Hyp1,Hyp2,flag_Hyp = 0,seed,standardize)
  #BFj0(BFji)
  bfj0_Hj<-bfju_Hj/bf0u_Hj
  }
  bfj0_Hj<-na.omit(bfj0_Hj)
  T0<-length(bfj0_Hj)
  # medianbfj0_Hj<-median(bfj0_Hj)
  kk<-lapply(1:T0, function (x) bfj0_Hj[x]>BFthresh )
  eta2<-sum(unlist(kk))/T0
   # kk<-lapply(1:T0, function (x) bfj0_Hj[x]<1/3 )
  # M_Hj<-sum(unlist(kk))/T0
  # kk<-lapply(1:T0, function (x) bfj0_Hj[x]>1/3&bfj0_Hj[x]<3 )
  # W_Hj<-sum(unlist(kk))/T0
  # CI_Hj<-round(quantile(bfj0_Hj,c(0.05,0.1,0.2,0.8,0.9,0.95)),2)
  #cat(paste('Median BFj0 = ',as.character(medianbfj0_Hj),'\n',sep=''))

  # kk<-lapply(1:T0, function (x) bfj0_Hj[x]>10 )
  # E_Hj<-sum(unlist(kk))/T0
  return(list(eta1,eta2))
}
