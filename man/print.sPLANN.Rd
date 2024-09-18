\name{print.sPLANN}
\alias{print.sPLANN}

\title{
Print Method for a Neural Network Based on the PLANN Method
}

\description{
This function prints a description of the neural network using the PLANN method. 
}

\usage{
\method{print}{sPLANN}(x, ...)
}

\arguments{
\item{x}{The result of a sPLANN fit.}
\item{...}{Further arguments passed to or from other methods.}
}

\value{No return value for this S3 method.}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}

\examples{
data(dataK)

splann <- sPLANN(Surv(time, event) ~ sex + stade + delay, data=dataK, inter=365.241, 
                          size=32, decay=0.01, maxit=100, MaxNWts=10000)

print(splann)
}


\keyword{Printing}
