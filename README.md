survivalPLANN: an R Package for SurvivalNeural Network
================

## Description

The R package ‘survivalPLANN’ contains a variety of functions to estimate a survival
neural network. The Partial Logistic Artificial Neural Networks (PLANN) are
implemented, proposed by Biganzoli et al. (1998). S3 methods are included to evaluate the
predictive capacities, as well as predictiions from new observations.

## Basic Usage

``` r
data(dataK)

splann <- survivalPLANN(Surv(time, event) ~ sex + stade + delay, data=dataK, inter=30, 
                          size=32, decay=0.01, maxit=200, MaxNWts=10000)

dnew <- data.frame(sex=c(1,2), delay=c(0,0), stade=c(0,0))

pred <- predict(splann, newdata = dnew)

# Predictions for a men or a women with no delay at the diagnostic of non-agressive cancer

plot(c(0,pred$times/365.241), c(1,pred$predictions[1,]), ylab="Patient survival",
  xlab="Post-diagnosis time in years", type="l")
  
lines(c(0,pred$times/365.241), c(1,pred$predictions[2,]), type="l", col=2)
```


## Installation

To install the latest release from CRAN:

``` r
install.packages("survivalPLANN")
```

To install the development version from GitHub:

``` r
remotes::install_github("chupverse/survivalPLANN")
```

## Reporting bugs

You can report any issues at this
[link](https://github.com/chupverse/survivalPLANN/issues).
