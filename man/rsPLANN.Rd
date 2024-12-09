\name{rsPLANN}
\alias{rsPLANN}

\title{
Relative Survival Neural Network by Using the PLANN Method
}
\description{
This function allows to compute a predicting tool for relative survival additive risks with neural network using the PLANN method.
}

\usage{
rsPLANN(formula, data, pro.time=NULL, inter, size = 32, decay = 0.01,
                    maxit = 100, MaxNWts = 10000, trace = FALSE,
                    ratetable)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the
right. The response must be a survival object as returned by the \code{Surv} function. On the right
side, it is necessary to specify a '+\code{ratetable}(age = age, year = year, sex = sex)' term. The
time must be in days.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{pro.time}{A numeric value with the maximum prognostic time in days. If \code{NULL}, the maximum prognostic time is the last time for wich an event is observed.}
\item{inter}{A numeric value representing the length of the intervals.}
\item{size}{A numeric value for the number of units in the hidden layer. Default is set to 32}
\item{decay}{A numeric value for the parameter for weight decay. Default is set to 0.01}
\item{maxit}{A numeric value for the maximum number of iterations. Default is set to 100.}
\item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the code, but
increasing MaxNWts will probably allow fits that are very slow and time-consuming. Default is set to 10000}
\item{trace}{A logical argument for tracing optimization. Default is set to FALSE.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, year and sex.}
}

\value{
\item{formula}{The \code{formula} object.}
\item{data}{The data frame given with the computed \code{Intervals} column added.}
\item{ratetable}{The ratetable object given to the function.}
\item{ays}{A data.frame of the age, year and sex variables given in the \code{ratetable} term.}
\item{pro.time}{The maximum prognostic time.}
\item{fitsurvivalnet}{The fitted \code{nnet} object.}
\item{times}{A vector with the \code{P} prognostic times up to \code{pro.time}}
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

\details{
The Biganzoli's PLANN method is used to predict the overall survival.
}

\references{Biganzoli E, Boracchi P, Mariani L, and et al. Feed forward neural networks for the 
analysis of censored survival data: a partial logistic regression approach. Stat Med, 17:1169-86,
1998.}


\examples{
data(dataK)
data(fr.ratetable)

rs1 <- rsPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker + ratetable(age=age,                year = year, sex=sexchara), data = dataK[1:200,],
               pro.time = 365.241*6, inter=365.241, size = 16, decay = 0.01,
               maxit = 100, MaxNWts=10000, ratetable=fr.ratetable
               )

plot(rs1$times/365.241, rs1$mpredictions$CIF_C,
xlab="Post-diagnostic time (years)",
     ylab="Cumulative incidence functions",
     ylim=c(0,1), col="red", type="s",)

lines(rs1$times/365.241, rs1$mpredictions$CIF_P, type="s", col="blue")

legend("topleft", legend=c("Cancer-related death", "Death due to other causes"), col=c("red", 
"blue"), lty=c(1,1))
}

\keyword{Neural Network}
