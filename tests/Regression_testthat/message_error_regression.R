#test-that the error message in ANOVA

library(testthat)
#========================================================================================")
#  Test the error messages if BFthresh is smaller than 1 or eta is smaller than 0
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='Ha',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=0.5,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='Ha',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=-0.5,seed=10,standardize=FALSE,ratio=c(1,1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='beta1>0&beta2>0&beta3>0',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,R_square2=0.13,T_sim=10000,
                                                     BFthresh=0.3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='beta1>0&beta2>0&beta3>0',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,R_square2=0.13,T_sim=10000,
                                                     BFthresh=3,eta=-0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})

#========================================================================================")
#  Test the error messages if the Hypothesis or the ratio is wrong
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1<beta2=beta3=0',Hyp2='Ha',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3',Hyp2='beta1>beta2>beta3',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),
                                                     R_square1=0,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=TRUE,ratio=c(-1,2,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2<beta3=0',Hyp2='Ha',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})

test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1>0&beta2<0&beta3>0',Hyp2='Hc',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1>0&beta2>0&beta3>0',Hyp2='Hc',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,-1,1)))})
test_that("SSDbain",{expect_error(Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='Ha',k=3,
                                                     rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0.13,
                                                     R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1)))})

