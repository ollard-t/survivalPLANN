\name{plot.sPLANN}
\alias{plot.sPLANN}

\title{
Calibration Plot for a Survival Neural Network Using the PLANN Method.
}
\description{
A calibration plot of an object of the class \code{sPLANN}.
}

\usage{
\method{plot}{sPLANN}(x, n.groups=5, pro.time=NULL,
newdata=NULL, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{sPLANN}.}
   \item{n.groups}{A numeric value with the number of groups by their class probabilities. The 
   default is 5.}
  \item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}: the calibration plot is performed from the same subjects of the training sample.}
  \item{...}{Additional arguments affecting the plot.}
}

\details{
The plot represents the survival obtained by the Kaplan-Meier estimator and the related 95\% 
confidence intervals, against the mean of the predictive values for individuals stratified into 
groups of the same size according to the percentiles. The identity line is usually included for reference.
}

\value{No return value for this S3 method.}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The multivariate PH regression with a Weibull distribution

model <- sPLANN(Surv(time, event) ~ biomarker + sex + stade + delay 
                , data = dataK, inter = 365.241)

# the calibration plots for predictions at 2, 4, 6 and 8 years

par(mfrow=c(2, 2))

plot(model, n.groups=2, pro.time=2*365.24,
  main="Predction at 2 years")

plot(model, n.groups=2, pro.time=4*365.24, 
  main="Predction at 4 years")

plot(model, n.groups=2, pro.time=6*365.24, 
  main="Predction at 6 years")

plot(model, n.groups=2, pro.time=8*365.24, 
  main="Predction at 8 years")
}

\seealso{
\code{\link{plot.default}}, \code{\link{sPLANN}}.
}

\keyword{Plotting}
