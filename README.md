# SSDbain
SSDbain was built under R version 3.6.3.    
The function SSDttest is used to compute the sample size required per group for the Bayesian t-test and Bayesian Welch's test.
The function SSDANOVA and SSDANOVA_robust in the R package SSDbain computes the sample size for the Bayesian ANOVA, Welch's ANOVA, and robust ANOVA.
The function SSDRegression is used to compute the sample size required per group for the linear regression models.

## Installation
Install the latest release version of `SSDbain` from github:

### install devtools
```
install.packages("devtools")
```

### Load devtools package for install_github()
```
library(devtools)

```
### get SSDbain from github
```
install_github("Qianrao-Fu/SSDbain",upgrade="never")

```
### Load SSDbain package
```
library(SSDbain)
```

## Citing SSDbain

You can cite the R-package with the following citation:

 Q.Fu. (2021). SSDbain: Sample Size Determination using AAFBF for Bayesian Hypothesis Testing Implemented in bain. (Version 0.1.0), R package.

## Contributing and Contact Information

If you have ideas, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.
Contributions in code must adhere to the [tidyverse style
guide](https://style.tidyverse.org/).

-   File a GitHub issue [here](https://github.com/Qianrao-Fu/SSDbain)
-   Make a pull request [here](https://github.com/Qianrao-Fu/SSDbain/pulls)

By participating in this project, you agree to abide by the [Contributor
Code of Conduct v2.0](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html).
