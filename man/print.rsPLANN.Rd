\name{print.rsPLANN}
\alias{print.rsPLANN}

\title{
Print Method for a Neural Network Based on the PLANN Method for relative survival
}

\description{
This function prints a description of the neural network using the PLANN method for relaive survival. 
}

\usage{
\method{print}{rsPLANN}(x, ...)
}

\arguments{
\item{x}{The result of a rsPLANN fit.}
\item{...}{Further arguments passed to or from other methods.}
}

\value{No return value for this S3 method.}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}

\examples{
data(dataK)
data(fr.ratetable)

rsplann <- rsPLANN(Surv(time, event) ~ sex + stade + delay + ratetable(age =age,      
           year = year, sex = sexchara), data=dataK, inter=365.241,
          size=32, decay=0.01, maxit=100, MaxNWts=10000, ratetable=fr.ratetable)

print(rsplann)
}


\keyword{Printing}
