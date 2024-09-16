\name{predict.sPLANN}
\alias{predict.sPLANN}

\title{
Predict Survival From a Neural Network Based on the PLANN Method
}

\description{
This function produces survival prediction from a neural network based on the PLANN method. 
}

\usage{
\method{predict}{sPLANN}(object, newdata = NULL, newtimes = NULL, ...)
}

\arguments{
\item{object}{The result of the \code{sPLANN} function.}
\item{newdata}{An optional data frame comprising of new examples to be predicted. If \code{NULL}, the data frame used is the one used for training in the \code{sPLANN} function.}
\item{newtimes}{A optional numeric vector comprising of times to get survival estimations. If \code{NULL}, the times are the intervals used in the \code{sPLANN} function.}
\item{...}{Further arguments passed to or from other methods.}
}


\value{
\item{times}{The times used for the predicitions.}
\item{predictions}{A data frame comprising of the survival predictions from the neural network.}
}


\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}

\examples{
data(dataK)

splann <- sPLANN(Surv(time, event) ~ sex + stade + delay, data=dataK, inter=365, 
                          size=32, decay=0.01, maxit=200, MaxNWts=10000)

dnew <- data.frame(sex=c(1,2), delay=c(0,0), stade=c(0,0))

pred <- predict(splann, newdata = dnew)

# Predictions for a men or a women with no delay at the diagnostic of non-agressive cancer

plot(c(0,pred$times/365.241), c(1,pred$predictions[1,]), ylab="Patient survival",
  xlab="Post-diagnosis time in years", type="l")
lines(c(0,pred$times/365.241), c(1,pred$predictions[2,]), type="l", col=2)
}

\keyword{Predicting}
