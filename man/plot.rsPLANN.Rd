\name{plot.rsPLANN}
\alias{plot.rsPLANN}

\title{
Calibration Plot for a Relative Survival Neural Network Using the PLANN Method.
}
\description{
A calibration plot of an object of the class \code{rsPLANN}.
}

\usage{
\method{plot}{rsPLANN}(x, n.groups=5, type = "relative", pro.time=NULL,
newdata=NULL, ...)
}

\arguments{
   \item{x}{An object returned by the function \code{rsPLANN}.}
   \item{n.groups}{A numeric value with the number of groups by their class probabilities. The 
   default is 5.}
   \item{type}{The type of predictions to be used in the calibration plot. Can be set to "relative",
   "net", "CIF_C" or "CIF_P" (see \code{details}). By default, \code{type} is set to "relative".}
   \item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities.}
   \item{newdata}{An optional data frame containing the new sample for validation with covariate 
   values, follow-up times, and event status. The default value is \code{NULL}: the calibration plot
   is performed from the same subjects of the training sample.}
   \item{...}{Additional arguments affecting the plot.}
}

\details{
The plot represents the Relative or Net Survival predictions or the Cumulative Incidence function 
of the cause specific or population obtained by different estimators and the related 95\% 
confidence intervals, against the mean of the predictive values for individuals stratified into 
groups of the same size according to the percentiles. The identity line is usually included for 
reference. If \code{type} is set to "relative", the relative survival predictions of the rsPLANN 
object are compared to the Ederer estimator. If set to "net", the net survival predictions are
compared to the Pohar-Perme estimator. The Ederer and the Pohar-Perme estimations are from the 
\code{rs.surv} function from the \code{relsurv} package. If set to "CIF_C" or "CIF_P", the cause 
specific or population Cumulative Incidence Function are compared to the ones computed by the 
\code{cmp.rel} function from the \code{relsruv} package. 
}

\value{No return value for this S3 method.}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The multivariate PH regression with a Weibull distribution

model <- rsPLANN(Surv(time, event) ~ biomarker + sex + stade + delay + ratetable(age=age ,year = 
year, sex=sexchara), data = dataK, inter = 365.241, ratetable = fr.ratetable)

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
\code{\link{plot.default}}, \code{\link{rsPLANN}}.
}

\keyword{Plotting}
