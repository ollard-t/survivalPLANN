\name{survivalPLANN}
\alias{survivalPLANN}

\title{
Survival Neural Network by Using the PLANN Method
}
\description{
This function allows to compute a neural network using the PLANN method.
}

\usage{
survivalPLANN(formula, data, inter, size = 32, decay = 0.01,
                          maxit = 100, MaxNWts = 10000, trace = FALSE ...)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{inter}{A numeric value representing the length of the intervals.}
\item{size}{A numeric value for the number of units in the hidden layer. Default is set to 32}
\item{decay}{A numeric value for the parameter for weight decay. Default is set to 0.01}
\item{maxit}{A numeric value for the maximum number of iterations. Default is set to 100.}
\item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the code, but increasing MaxNWts will probably allow fits that are very slow and time-consuming. Default is set to 10000}
\item{trace}{A logical argument for tracing optimization. Default is set to FALSE.}
\item{...}{Further arguments passed to or from other methods.}
}

\value{
\item{formula}{The \code{formula} object.}
\item{fitsurvivalnet}{The fitted \code{nnet} object.}
\item{data}{The data frame given with the computed \code{Intervals} column added.}
\item{data_dup}{The base data frame with the duplicated rows based on the \code{Intervals} column. }
\item{call}{The function call for the creation of the \code{nnet} object.}
\item{inter}{The interval length.}
\item{size}{The number of units in the hidden layer.}
\item{decay}{The value of the parameter for weight decay.}
\item{maxit}{The value of the parameter for maxit.}
\item{MaxNWts}{The value of the parameter for MaxNWts.}
\item{coefnames}{The name(s) of the covariate(s) used in the model.}
\item{y}{The values of \code{time} and \code{event}.}
\item{x}{The model covariates values.}
\item{intervals}{The intervals limits values.}
\item{missing}{The rows with a missing value in the data base.}
}

\details{
The \code{nnet} loss function parameter is fixed to \code{entropy} as Biganzoli's PLANN method minimize the cross-entrophy loss function which correspondS to maximizing the loglikelihood.
}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}


\examples{
data(dataK)

splann <- survivalPLANN(Surv(time, event) ~ sex + stade + delay, data=dataK, inter=365, 
                          size=32, decay=0.01, maxit=100, MaxNWts=10000)

print(splann)
}

\keyword{Neural Network}
