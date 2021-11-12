beta_Hc2<-function(k,beta){
  library(gtools)
  x <- c(1:k)
  #pick 2 balls from the urn with replacement
  #get all permutations
  database=permutations(n=k,r=k,v=x)
  num=nrow(database)

  database_temp<-list()
  database_temp[[1]]=database[1,]
  num_total<-sum(1:k-1)+1
  j<-1
  while(j<num_total){
    num_j<-length(database_temp[[j]])/k
    matrix_temp<-matrix(database_temp[[j]],num_j,k)
    vector_diff<-c()
    x<-1
    flag_row<-c()

    for (i in 1:num_j)
    {
      for (ii in 1:num)
      {
        vector_diff<-database[ii,]-matrix_temp[i,]
        index<-which(vector_diff!=0)
        if(length(index)==2)
        {
          if((index[2]-index[1])==1)
          {flag_row[x]<-ii
          x<-x+1}
        }
      }
    }
    database_temp[[j+1]]<-unique(database[flag_row,])
    j<-j+1
  }

  data_final<-rbind(database_temp[[1]],database_temp[[2]])
  if(k>2){
  for (i in 3:num_total)
  {
    data_final<-rbind(data_final,database_temp[[i]])
  }}
  data_final<-unique(data_final)
  seq<-data_final[median(2:num),]
  return(beta[seq])}

beta_Hc1<-function(k,beta){
  Num_mid<-2^(k-1)
  i<-1
  sum<-choose(k, i)
  while(sum<Num_mid){
    i<-i+1
    sum<-sum+choose(k,i)
  }
  j<-i
  beta_Hc<-beta
  for (i in 1:j)
  {
    beta_Hc[i]=-beta[i]
  }
  return(beta_Hc)}

cal_beta<-function(k,R_square1,R_square2,Hyp1,Hyp2,ratio,rho){
  varnames<-c()
  b<-character()
  for (i in 1:k){
    varnames[i]=paste('beta',as.character(i),sep='')
  }

  num_eq <- sum(unlist(gregexpr("[=]",Hyp1))>0)
  num_0 <- sum(unlist(gregexpr("[0]",Hyp1))>0)
  num_larger <- sum(unlist(gregexpr("[>]",Hyp1))>0)
  if(num_eq==k&&num_0>0)
    beta1<-numeric(k)
  else if(num_eq==k-1&&num_0==0){
    beta<-numeric(k)+1
    beta_r<-matrix(runif(k*k),nrow = k)
    for (i in 1:k){
      for (j in 1:k){
        beta_r[i,j]<-beta[i]*beta[j]*rho[i,j]
      }
    }
    d<-sqrt(R_square1/sum(beta_r))
    beta1<-d*beta
  }
  else{
    if(R_square1==0)
      beta1<-rep(0,k)
    else
      beta1<-cal_beta_general(k,R_square1,Hyp1,ratio,rho)
  }

  if(Hyp2=='Hc')
  {
     if(num_larger==k&&num_0==k){
       beta2<-beta_Hc1(k,beta1)}
     else if(num_larger==k-1){
       beta2<-beta_Hc2(k,beta1)
     }
     else
      #cat(paste('Please check your input hypothesis!  ','\n',sep=''))
      beta2<-numeric(k)
  }
  else if(Hyp2=='Ha'){
    beta<-numeric(k)
    beta<-ratio
   # beta=seq(from=k,to=1,by=-1)
    beta_r<-matrix(runif(k*k),nrow = k)
    for (i in 1:k){
      for (j in 1:k){
        beta_r[i,j]<-beta[i]*beta[j]*rho[i,j]
      }
    }
    d<-sqrt(R_square2/sum(beta_r))
    beta<-d*beta
    beta2<-beta
  }else{
    num_eq <- length(unlist(gregexpr("[=]",Hyp2)))
    num_0 <- length(unlist(gregexpr("[0]",Hyp2)))
    num_larger <- length(unlist(gregexpr("[>]",Hyp2)))
    if(num_eq==k&&num_0>0)
      beta2<-numeric(k)
    else
      beta2<-cal_beta_general(k,R_square2,Hyp2,ratio,rho)
  }
  return(list(beta1,beta2))
}

cal_beta_general<-function(k,R_square,Hyp,ratio,rho){
  varnames<-c()
  for (i in 1:k){
    varnames[i]=paste('beta',as.character(i),sep='')
  }

  # beta<-c(1:k)*0
  # EI_matrix1<-matrix_trans(varnames,Hyp)
  # ERr1<-EI_matrix1[[1]]
  # IRr1<-EI_matrix1[[2]]
   matches <- regmatches(Hyp, gregexpr("[[:digit:]]+", Hyp))
  index<-as.numeric(unlist(matches))
  # if(length(which(index==0))==0){
  #   beta[index]=seq(from=k,to=1,by=-1)
  #   }else if(length(which(index==0))==1){
  # NE<-index[(which(index==0)+1):(k+1)]
  # beta[NE]=seq(from=-1,to=-length(NE),by=-1)
  # beta[index[1:which(index==0)-1]]=seq(from=which(index==0)-1,to=1,by=-1)}else{}
  #
  # if(nrow(IRr1)>=k&&ncol(IRr1)>=k){
  #   if(sum(IRr1[1:k,1:k]-diag(x=1,k,k))==0){
  #     beta=seq(from=5*k-4,to=1,by=-5)
  #   }}
  #
  beta<-ratio
  beta_r<-matrix(runif(k*k),nrow = k)
  for (i in 1:k){
    for (j in 1:k){
      beta_r[i,j]<-beta[i]*beta[j]*rho[i,j]
    }
  }
  d<-sqrt(R_square/sum(beta_r))
  beta[index]<-d*beta
  return(beta)
}



matrix_trans<-function(varnames,hyp){
  Rr1<-create_matrices(varnames, hyp)
  Rr1$ERr1[is.null(Rr1$ERr1)] <- 0
  nrow_ERr1<-nrow(Rr1$ERr1)
  ncol_ERr1<-ncol(Rr1$ERr1)
  nrow_ERr1[is.null(nrow_ERr1)] <- 0
  ncol_ERr1[is.null(ncol_ERr1)] <- 0

  Rr1$IRr1[is.null(Rr1$IRr1)] <- 0
  nrow_IRr1<-nrow(Rr1$IRr1)
  ncol_IRr1<-ncol(Rr1$IRr1)
  nrow_IRr1[is.null(nrow_IRr1)] <- 0
  ncol_IRr1[is.null(ncol_IRr1)] <- 0


  ERr1<-matrix(Rr1$ERr1,nrow_ERr1,ncol_ERr1)
  IRr1<-matrix(Rr1$IRr1,nrow_IRr1,ncol_IRr1)
  return(list(ERr1,IRr1))
}
