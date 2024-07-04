\name{predict.survivalPLANN}
\alias{predict.survivalPLANN}

\title{
Predict Method for survivalPLANN objects
}
\description{
This function produce a prediction from a fitted neural network using the PLANN method to estimate survival. 

}
\usage{
\method{predict}{survivalPLANN}(object, newdata = NULL, newtimes = NULL, ...)
}

\arguments{
\item{object}{The result of a survivalPLANN fit.}
\item{newdata}{A data frame comprising of new examples to be predicted. If \code{NULL}, the data frame used is the one given in the \code{survivalPLANN} object.}
\item{newtimes}{A numeric vector comprising of times to get survival estimations. If \code{NULL}, the time values chosen are the ones from the \code{survivalPLANN} object.}
\item{...}{Further arguments passed to or from other methods.}

}


\value{
\item{times}{The times used for the predicitions.}
\item{predictions}{A data frame comprising of the predictions from the neural network.}

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
                          
db.new <- data.frame(X1=c(0,0,1,1), X2=c(0,1,0,1))

predict(splann, newdata = db.new)
}


\keyword{Neural Network}
