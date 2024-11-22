\name{predict.rsPLANN}
\alias{predict.rsPLANN}

\title{
Predictions From a Relative Neural Network Based on the PLANN Method
}

\description{
This function produces several predictions from a relative neural network based on the PLANN method. 
}

\usage{
\method{predict}{rsPLANN}(object, newdata = NULL, newtimes = NULL, ...)
}

\arguments{
\item{object}{The result of the \code{rsPLANN} function.}
\item{newdata}{An optional data frame comprising of new examples to be predicted. If \code{NULL}, the data frame used is the one used for training in the \code{rsPLANN} function.}
\item{newtimes}{A optional numeric vector comprising of times to get survival estimations. If \code{NULL}, the times are the intervals used in the \code{rsPLANN} function.}
\item{...}{Further arguments passed to or from other methods.}
}


\value{
\item{times}{The \code{P} times used for the predicitions.}
\item{ipredictions}{A list with the predicted outcomes for the \code{N} individals in \code{data}:
\code{survival_O} is the observed survival, \code{survival_P} is the population survival, i.e the 
survival of the individuals for individuals with the same caracteristics in the ratetable, 
\code{survivalR} is the relative survival and \code{survivalE} is the net survival. \code{CIF_C} is 
a \code{NxP} matrix with the cumulative probabilities of death caused by the disease, \code{CIF_P}
is a \code{NxP} matrix with the cumulative probabilities of death of the reference population in 
\code{ratetable}, \code{maxCIF_P} is a vector with (asymptotic) probabilty of being cured, 
\code{cure} is a \code{NxP} matrix with the probabilities of being cured conditionnaly to be alive 
at the progostic times.}
\item{mpredictions}{A list with the mean predicted outcomes : 
\code{survival_O} is the mean observed survival, \code{survival_P} is the mean population survival, 
i.e the survival of the individuals for individuals with the same caracteristics in the ratetable, 
\code{survivalR} is the mean relative survival and \code{survivalE} is the mean net survival. 
\code{CIF_C} is a vector with the mean cumulative probabilities of death caused by the disease,
\code{CIF_P} is a vector with the mean cumulative probabilities of death of the reference population
in \code{ratetable}. \code{cure} is a vector with the mean probabilities of being cured 
conditionnaly to be alive at the progostic times. \code{hazard_P} is the mean instaneous hazard in 
the population and \code{hazard_C} is the mean instaneous hazard in excess.}
}


\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86, 1998.}

\examples{
data(dataK)
data(fr.ratetable)

rs1 <- rsPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker + ratetable(age=age, year = 
year, sex=sexchara), data = dataK[1:200,],
               pro.time = 365.241*6, inter=365.241, size = 16, decay = 0.01,
               maxit = 100, MaxNWts=10000, ratetable=fr.ratetable)


plot(rs1$times/365.241, rs1$mpredictions$CIF_C,
xlab="Post-diagnostic time (years)",
     ylab="Cumulative incidence functions",
     ylim=c(0,1), col="red", type="s",)

lines(rs1$times/365.241, rs1$mpredictions$CIF_P, type="s", col="blue")

legend("topleft", legend=c("Cancer-related death", "Death due to other causes"), col=c("red", 
"blue"), lty=c(1,1))


p <- predict(rs1, newdata=dataK[1:200,], newtimes=365.241*c(3,5.5))
   
points(p$times/365.241, p$mpredictions$CIF_C, col="green", pch=3)
}

\keyword{Predicting}
