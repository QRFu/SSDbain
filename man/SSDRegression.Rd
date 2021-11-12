\name{SSDRegression}
\alias{SSDRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample size determination for the Bayesian Multiple Linear Regression with bain}
\description{
The function \code{SSDRegression} in the \code{R} package \code{SSDbain} computes the sample size for multiple linear regression using Bayes factors implemented in \code{bain} with informative hypothesis. \code{SSDbain} can be found at \url{https://github.com/Qianrao-Fu/SSDbain}. In \href{https://doi.org/10.31234/osf.io/3tr5f}{Fu (unpublished)}, it is elaborated how to use \code{SSDRegression}. Users are well advised to read this paper before using \code{SSDRegression}. Users are required to install package \code{bain} from CRAN before using \code{SSDbain} package.}

\usage{
SSDRegression(Hyp1,Hyp2,k,rho,R_square1,R_square2,T_sim = 1000,BFthresh,eta,seed = 10,standardize = TRUE, ratio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{\code{Hyp1}}{strings that specify the former hypothesis of a pair of interested hypotheses. For example, if \code{H0:"beta1=beta2=beta3=0"} vs \code{H1:"beta1>0&beta2>0&beta3>0"}, then \code{Hyp1:"beta1=beta2=beta3=0"}.}

\item{\code{Hyp2}}{strings that specify the latter hypothesis of a pair of interested hypotheses. For example, as the previous example, then \code{Hyp2:"beta1>0&beta2>0&beta3>0"}. Note: when Hyp2 is an unconstrained hypothesis, then \code{Hyp2="Ha"}; when hyp2 is the complement hypothesis of Hyp1, then \code{Hyp2="Hc"}.}

\item{\code{k}}{a positive integer that specifies the number of predictors in the multiple linear regression. It should be noted that the number of predictors in the hypotheses is equal to the number of predictors in the linear regression. For example, if the model is \eqn{yi=beta1 xi1+beta2 xi2+beta3 xi3+\epsilon i}, \code{k=3}.}

\item{\code{rho}}{a symmetric matrix that specifies the correlation between each pair of predictors. The diagonal of this matrix are ones.}

\item{\code{R_square1}}{parameter used to specify the coefficient of determination under \code{hyp1}.}

\item{\code{R_square2}}{parameter used to specify the coefficient of determination under \code{hyp2}.}

\item{\code{BFthresh}}{a numeric value not less than 1 that specifies the required size of the Bayes factor for the true hypothesis. For example, if \code{BFthresh=3}, and H0 is compared with H1, BF01\eqn{\ge}3 and BF10\eqn{\ge}3 should be satisfied.}

\item{\code{eta}}{a numeric value that specifies the probability that the Bayes factor is larger than \code{BFthresh} if either of the competing hypotheses is true.
}

\item{\code{T_sim}}{a positive integer that specifies the number of data sets sampled from the populations corresponding to the two competing hypotheses. \bold{It should be noted that at least \code{T_sim=10000} is required to guarantee the stability of the results.}}

\item{\code{seed}}{a positive integer that specifies the seed of \code{R}'s random numer of generator.
}

\item{\code{standardize}}{a logical value that specifies whether hypotheses regarding standardized or unstandardized regression coefficients are evaluated. With \code{standardize=TRUE}hypothese with respect to standardized regression coefficients are evaluated. With \code{standardize=FALSE} hypothese with respect to unstandardized regression coefficients are evaluated.
}

\item{\code{ratio}}{an optional vector that specifies the ratio between each pair of regression coefficients for each population in which the corresponding hypothesis is true.}
}

\note{The following competing hypothese can be executed:

1) \code{Hyp1='beta1=beta2=...=betaK=0'} vs \code{Hyp2='Ha'}, where \code{Ha} is the unconstrianed hypothesis. The default ratio of the regression coefficienst for \code{Ha} is 1,1,...,1.

2) \code{Hyp1='beta1=beta2=...=betaK=0'} vs \code{Hyp2='beta1>0,beta2>0,...,betaK>0'}, where each of the regression coefficient may also be smaller than zero. That is, both "<" and ">" can exist in a hypothesis. The default ratio of the regression coefficienst for \code{Hyp2} is 1,1,...,1.

3) \code{Hyp1='beta1=beta2=...=betaK'} vs \code{Hyp2='beta1*>beta2*>...>betaK*'}, where 1*, 2*, ..., K* are a re-ordering of the numbers 1, 2, ..., K. The ratio of the
regression coefficients 1,1, ... ,1 is used for \code{Hyp1} and the regression coefficients are computed such that R^2 = 0.13, and the ratio of the regression coefficients Kd, (K-1)d, ..., 3d, 2d, d is used for \code{Hyp2}.

4) \code{Hyp1='beta1>0, beta2>0, ..., betaK>0'}  vs \code{Hyp2='Hc'}, where Hc is the complement hypothesis of \code{Hyp1}. The ratio for \code{Hyp1} is consistent with \code{Hyp2} in 2), and the ratio is reordered using the representative hypothesis (see Appendix B in Fu, unpublished) for \code{Hyp2}. It should be noted that in this situation only ">" or "<" is allowed.

5) \code{Hyp1='beta1*>beta2*>...>betaK*'}  vs \code{Hyp2='Hc'}, where Hc is the complement hypothesis of \code{Hyp1}. The ratio for \code{Hyp1} is consistent with \code{Hyp2} in 3), and the ratio is reordered using the representative hypothesis (see Appendix B in Fu, unpublished) for \code{Hyp2}.

The standardized regression coefficients are used in 3) and 5) to ensure the regression coefficients are comparable.
}

\value{
The output resulting from analyses with \code{SSDRegression} contains:

1) required sample size N

2) the probability P(BFsv>BFthresh|Hs)

3) the probability P(BFvs>BFthresh|Hv)

Note: To perform the sensitive analysis, the results are provided for three different fractions b, 2b, and 3b, where b corresponds to fraction = 1 in the call to \code{bain}.

}

\references{
Fu, Q.(unpublished). Sample size determination for Bayesian testing of informative hypotheses in linear regression.
Gu, X., Mulder, J., and Hoijtink, H. (2017). Approximated adjusted fractional Bayes factors: A general method for testing informative hypotheses. British Journal of Mathematical and Statistical Psychology. \emph{British Journal of Mathematical and Statistical Psychology}, 71(2), 229-261. doi:10.1111/bmsp.12110

Hoijtink, H., Gu, X., and Mulder, J. (2018). Bayesian evaluation of informative hypotheses for multiple populations. \emph{British Journal of Mathematical and Statistical Psychology}. doi:10.1111/bmsp.12145
}
\author{
Qianrao Fu
}

\examples{


#Example 1: using the coefficient of determination R^2=0.13 to determine the sample size when H0 versus Ha.


Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='Ha',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1))


#Example 2: using the coefficient of determination R^2=0.13 to determine the sample size when H0 versus H1.


Res<-SSDRegression(Hyp1='beta1=beta2=beta3=0',Hyp2='beta1>0&beta2>0&beta3>0',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1))


#Example 3: using the coefficient of determination R^2=0.13 to determine the sample size when H0 versus H2.


Res<-SSDRegression(Hyp1='beta1=beta2=beta3',Hyp2='beta1>beta2>beta3',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=TRUE,ratio=c(3,2,1))


#Example 4: using the coefficient of determination R^2=0.13 to determine the sample size when H1 versus H1c.


Res<-SSDRegression(Hyp1='beta1>0&beta2>0&beta3>0',Hyp2='Hc',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0.13,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=FALSE,ratio=c(1,1,1))


#Example 5: using the coefficient of determination R^2=0.13 to determine the sample size when H2 versus H2c.


Res<-SSDRegression(Hyp1='beta1>beta2>beta3',Hyp2='Hc',k=3,rho=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),nrow=3),R_square1=0.13,R_square2=0.13,T_sim=10000,BFthresh=3,eta=0.8,seed=10,standardize=TRUE,ratio=c(3,2,1))


}

