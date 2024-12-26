survivalPLANN: an R Package for Survival Neural Network
================

## Description

The R package ‘survivalPLANN’ contains a variety of functions to predictive survival rates from neural network. It also allows us to predictive relative survivals. The Partial Logistic Artificial Neural Networks (PLANN) are implemented, proposed by Biganzoli et al. (1998). S3 methods are included to evaluate the predictive capacities, as well as predictions from new observations.

## Basic usage to predict the overall survival

``` r
## import libraries

library(survivalPLANN)
library(relsurv)
#library(survivalNET)
#library(lubridate)


## data management

data(dataK) # import the data base (colorectal cancers from the 'relsurv' package)

dataK$agey <- dataK$age/365.241

## estimation of the hyperparameters

pro.time <- floor(max(dataK$time[dataK$event==1])/365.241)  # 13 years

tune.sPLANN <- cvPLANN(Surv(time, event)~ sex + agey + stade + delay, data=dataK, cv=10,
                 pro.time=pro.time*365.241, inter=365.241/12, size=2:5, decay=c(0.01, 0.1))

tune.sPLANN$optimal$size
# [1] 2

tune.sPLANN$optimal$decay
# [1] 0.01

## estimation of the network according to the previous optimal parameters

splann <- sPLANN(Surv(time, event)~ sex + agey + stade + delay, data=dataK, 
                     pro.time=pro.time*365.241, inter=365.241/12, size=2, decay=0.01, maxit=1000)

# predictions for a 50-years old patientwith no delay at the diagnostic 
# of a non-agressive cancer according to the gender

dnew <- data.frame(sex=c(1,2), agey=c(50,60), stade=c(0,0), delay=c(0,0))

datap <- predict(splann, newdata = dnew) #survival predictions for dnew

plot(c(0,datap$times/365.241), c(1,datap$predictions[1,]), ylab="Patient survival",
  xlab="Post-diagnosis time in years", type="l") # sex=1 (male)
  
lines(c(0,datap$times/365.241), c(1,datap$predictions[2,]), type="l", col=2) #sex=2 (female)
```


## Individual and marginal predictions of overall and relative survival rates

``` r
data("fr.ratetable") # import the table with the expected population mortality

datap <- predictRS(object=splann, data=dataK,
                 ratetable=fr.ratetable, age="age", sex="sexchara", year="year")

# the predicted overall survival curves of the first 100 patients

plot(survfit(Surv(time/365.241, event) ~ 1, data = dataK),
     ylab="Overall survival", xlab="Time (years)", conf.int = FALSE,
     lwd=2, col="red")
     
for (i in 1:100) {
lines(x=datap$times/365.241, y=datap$ipredictions$overall_survival[i,],
     col="gray", type="s") }

legend("topright", c("Kaplan-Meier estimator", "Individual predictions"),
    col=c("red", "gray"), lty=c(1,1), lwd=c(2,1))
```


## Describing the performences of 'splann' to predict the overall survival

``` r
plot(survfit(Surv(time/365.241, event) ~ 1, data = dataK),
     ylab="Overall survival", xlab="Time (years)", 
     lwd=1, col="black") # the non-parametric Kaplan-Meier estimator
     
lines(x=datap$times/365.241, y=datap$mpredictions$overall_survival,
      col="red", type="s")

legend("topright", c("Kaplan-Meier estimator", "Mean of the individual predictions"),
       col=c("black", "red"), lty=c(1,1), lwd=c(1,1))
```


## Describing the performences of 'splann' to predict the net survival

``` r
fit_net <- rs.surv(Surv(time, event) ~ 1, data=dataK, ratetable=fr.ratetable,
                    rmap=list(age=age, sex= sex, year=year),
                    method = "pohar-perme") # the non-parametric Pohar-Perme estimator

plot(fit_net, col=1, lwd=1, lty=1, xscale = 365.241, xlab="Time (years)", ylim=c(0,1))

lines(x=datap$times, y=datap$mpredictions$net_survival,
      type="s", col="red")

legend("topright", c("Pohar-Perme estimator", "Mean of the individual predictions"),
    col=c("black", "red"), lty=c(1,1), lwd=c(1,1))
```


## Describing the performences of 'splann' to predict the cumulative incidence function

``` r
cmp_fit <- cmp.rel(Surv(time, event) ~ 1, data=dataK, ratetable=fr.ratetable,
                    rmap=list(age=age, sex= sex, year=year)) # the non-parametric Pohar-Perme estimator 

plot(cmp_fit, col=c("black", "black"), lty=c(1,2), xscale = 365.241,
     xlab="Time (years)", ylim=c(0,1))

lines(datap$times/365.241, datap$mpredictions$excess_cif, type="s", col="red", lty=1)

lines(datap$times/365.241, datap$mpredictions$population_cif, type="s", col="red", lty=2)
```


## Describing the performences of 'splann' to predict the relative survival ratio

``` r
fit_sr <- rs.surv(Surv(time, event) ~ 1, data=dataK, ratetable=fr.ratetable,
                  rmap=list(age=age, sex= sex, year=year), method = "ederer1")

plot(fit_sr, col=1, lwd=1, lty=1, xscale = 365.241, xlab="Time (years)", ylim=c(0,1))

lines(x=datap$times, y=datap$mpredictions$relative_ratio_survival,
      type="s", col="red")

legend("topright", c("Ederer-I estimator", "Mean of the individual predictions"),
    col=c("black", "red"), lty=c(1,1), lwd=c(1,1))
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
