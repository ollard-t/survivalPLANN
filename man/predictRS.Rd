\name{predictRS}
\alias{predictRS}

\title{
Predict Net Survivals From a Neural Network Based on the PLANN Method
}

\description{
This function produces several net survival predictions from a neural network based on the PLANN method. 
}

\usage{predictRS(object, data, ratetable, age, year, sex) }

\arguments{
\item{object}{The result of the \code{sPLANN} function.}
\item{data}{A data frame comprising of new observed predictors}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{A numeric vector with the patient age at baseline in days.}
\item{year}{A numeric vector with the baseline date in the study respecting the date format: number of days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
}


\value{  
\item{nnet}{The fitted survival neural network  from the \code{sPLANN} function.}
\item{time}{A vector with the post-baseline days for which the survival predictions are performed.}
\item{ipredictions}{A list with the individual predicted outcomes : COMPLETE THE LIST.}
\item{mpredictions}{A list with the marginal predicted outcomes : COMPLETE THE LIST.}
}


\references{Pohar Perme M, Stare J, Esteve J. On Estimation in Relative Survival. Biometrics, 68:113-120, 2012.}

\examples{
data(dataK) # import the data base (colorectal cancers from the 'relsurv' package)
dataK$agey <- dataK$age/365.241
pro.time <- floor(max(dataK$time[dataK$event==1])/365.241)  # 13 years

splann <- sPLANN(Surv(time, event)~ sex + agey + stade + delay, data=dataK, 
                     inter=365.241/12, size=2, decay=0.001, maxit=1000)

data("fr.ratetable") # import the table with the expected population mortality

# predictions of the relative survival for the first 10 patients
#datap <- predictRS(object=splann, data=dataK[1:10,],
#                 ratetable=fr.ratetable, age="age", sex="sexchara", year="year")
                 
#plot(x=datap$times/365.241, y=datap$ipredictions$relative_survival[1,],
#     col="gray", type="s", xlab="Time in years", ylab="Relative survival")
     
#for (i in 2:5) {
#     lines(x=datap$times/365.241, y=datap$ipredictions$relative_survival[i,],
#           col="gray", type="s")
#     }
}

\keyword{Predicting}
\keyword{Net survival}
\keyword{Relative survival}