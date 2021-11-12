#test-that the error message in ANOVA

library(testthat)
#========================================================================================")
#  Test the error messages if BFthresh is smaller than 1 or eta is smaller than 0
#========================================================================================")
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(1.5,0.75,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(1.5,0.75,0.75),BFthresh=0.5,eta=0.8,T=100,seed=10))})
test_that("SSDbain", {expect_error(SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(1.75,1.75,1.75),kurts=c(5.89,5.89,5.89),var=c(1.5,0.75,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain", {expect_error(SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(1.75,1.75,1.75),kurts=c(5.89,5.89,5.89),var=c(1.5,0.75,0.75),BFthresh=0.5,eta=0.8,T=100,seed=10))})


#========================================================================================")
#  Test the error messages if the mean values are not consistent with the hypothesis
#========================================================================================")
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="mu2>mu1>mu3",type="unequal",f1=c(0.2,0,3,0.5),f2=c(0.2,0,3,0.5),var=c(1.5,0.75,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="mu1>mu2>mu3",type="equal",f1=c(0.2,0,3,0.5),f2=c(0.2,0,3,0.5),var=c(1,1,1),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu2>mu1>mu3",hyp2="mu1>mu2>mu3",type="unequal",f1=c(0.2,0,3,0.5),f2=c(0.2,0,3,0.5),var=c(1.5,0.75,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu3>mu1>mu2",hyp2="mu1>mu2>mu3",type="equal",f1=c(0.2,0,3,0.5),f2=c(0.2,0,3,0.5),var=c(1,1,1),BFthresh=1,eta=-0.8,T=100,seed=10))})

#========================================================================================")
#  Test the error messages if the number of mean values are not consistent with the number of variables
#========================================================================================")
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=c(0.25,0.2),f2=c(0.25,0.3),var=c(1.5,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=c(0.25,0.25,0.25,0.2),f2=c(0.25,0.25,0.25,0.3),var=c(1.5,0.75,0.75,1),BFthresh=1,eta=-0.8,T=100,seed=10))})


#========================================================================================")
#  Test the error messages if the variance values are not consistent with its type
#========================================================================================")
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="equal",f1=0.25,f2=0.25,var=c(1.5,0.75,0.75),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(1,1,1),BFthresh=1,eta=-0.8,T=100,seed=10))})


#========================================================================================")
#  Test the error messages if the number of variance values are not consistent with the number of variables
#========================================================================================")
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(1.5,0.75,0.75,1),BFthresh=1,eta=-0.8,T=100,seed=10))})
test_that("SSDbain",{expect_error(SSDANOVA(hyp1="mu1>mu2>mu3",hyp2="Hc",type="unequal",f1=0.25,f2=0.25,var=c(0.75,1),BFthresh=1,eta=-0.8,T=100,seed=10))})


#========================================================================================")
#  Test the error messages if if kurtosis>skewness^2-2
#========================================================================================")
test_that("SSDbain", {expect_error(SSDANOVA_robust(hyp1="mu1=mu2=mu3",hyp2="mu1>mu2>mu3",f1=0,f2=0.25,skews=c(1.75,1.75,1.75),kurts=c(1,5.89,5.89),var=c(1.5,0.75,0.75),BFthresh=0.5,eta=0.8,T=100,seed=10))})
