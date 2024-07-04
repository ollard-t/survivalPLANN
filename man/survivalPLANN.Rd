\name{survivalPLANN}
\alias{survivalPLANN}

\title{
Neural Network survival estimation using the PLANN method
}
\description{
This function allows to compute a neural network using the PLANN method for 
modelling survival.
}
\usage{
survivalPLANN(formula, data, inter, size = 32, decay = 0.01,
                          maxit = 100, MaxNWts = 10000, ...)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{inter}{A numeric value representing the length of the intervals.}
\item{size}{A numeric value for the number of units in the hidden layer. Default is set to 32}
\item{decay}{A numeric value for the parameter for weight decay. Default is set to 0.01}
\item{maxit}{A numeric value for the maximum number of iterations. Default is set to 100.}
\item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the
code, but increasing MaxNWts will probably allow fits that are very slow and time-consuming. Default is set to 10000}
\item{...}{Further arguments passed to or from other methods.}

}


\value{ 
\item{formula}{The \code{formula} object.}
\item{fitsurvivalnet}{The fitted \code{nnet} object.}
\item{data}{The data frame given with the computed \code{Intervals} column added}
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
The \code{nnet} loss function parameter is fixed to \code{entropy} as Biganzoli's PLANN method minimize the cross-entrophy loss function which correspond to maximizing the loglikelihood.
}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}


\examples{
# database simulation

Fm1 <- function(u, x, coef, sigma, theta, nu)
{ sigma * ((1-log(u) *
              exp( as.matrix(x) \%*\% as.vector(-1*coef) ))^theta - 1)^(1/nu)
}

N <- 500; sigma <- 40; nu <- 2; theta <- 1; beta <- log(1.5)
X1 <- rbinom(N, 1, prob=0.50)
X2 <- rbinom(N, 1, prob=0.25)
X3 <- rnorm(N, 0, 1)
U <- runif(N, min=0, max=1)
t.sim <- as.vector(Fm1(u=U, x=cbind(X1, X2, X3), coef=c(beta, beta, beta),
                       sigma, theta, nu))
c.sim <- runif(n = N, min = 0, max = 60)
y.sim <- pmin(t.sim, c.sim)
db <- data.frame(event=1*(t.sim==y.sim), time=y.sim, X1, X2)

splann <- survivalPLANN(Surv(time, event) ~ X1 + X2, data=db, inter=2, 
                          size = 32, decay = 0.01, maxit =100, MaxNWts=10000)
}

\keyword{Neural Network}
