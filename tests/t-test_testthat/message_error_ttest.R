#test-that the error message in ANOVA

library(testthat)
#========================================================================================")
#  Test the error messages if BFthresh is smaller than 1 or eta is smaller than 0
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=0.5,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=0.5,eta=-0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=0.5,eta=0.8,Hypothesis='one-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=NULL,BFthresh=0.5,eta=-0.8,Hypothesis='one-sided',T=10000))})

#========================================================================================")
#  Test the error messages if the number of mean values is wrong
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.8,0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5),var=NULL,BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5),var=NULL,BFthresh=3,eta=0.8,Hypothesis='one-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.8,0.5,0),var=NULL,BFthresh=3,eta=0.8,Hypothesis='one-sided',T=10000))})


#========================================================================================")
#  Test the error messages if the variance values are not consistent with its type
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=c(1.3,0.8),BFthresh=0.5,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='unequal',Population_mean=c(0.5,0),var=c(1,1),BFthresh=0.5,eta=-0.8,Hypothesis='two-sided',T=10000))})


#========================================================================================")
#  Test the error messages if the number of variance is wrong
#========================================================================================")
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=c(1.5,0.7,0.7),BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=c(1.5),BFthresh=3,eta=0.8,Hypothesis='two-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=c(1.5,0.7,0.7),BFthresh=3,eta=0.8,Hypothesis='one-sided',T=10000))})
test_that("SSDbain",{expect_error(Res<-SSDttest(type='equal',Population_mean=c(0.5,0),var=c(1.5),BFthresh=3,eta=0.8,Hypothesis='one-sided',T=10000))})



