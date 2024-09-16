
rsPLANN <- function(formula, data, pro.time=NULL, inter, size= 32, decay=0.01,
                    maxit=100, MaxNWts=10000, trace=FALSE,
                    ratetable, age, year, sex) # @Thomas: inclure ratetable et ses arguments dans formula
{
  
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(inter)) stop("an inter argument is required")
  if (as.character(class(formula)) != "formula") stop("The formula argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (as.character(class(inter)) != "numeric") stop("The inter argument must be numeric")
  
  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
  ####### data management
  
  splann <- sPLANN(formula, data=data, pro.time=pro.time, inter=inter, 
                          size=size, decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  predO <- predict(splann, newtimes=splann$intervals)
  
  times <- predO$times
  
  exphaz <- function(x,age,sex,year) { expectedhaz(ratetable, age=age, sex=sex, year=year, time=x)}
  
  survO <- as.matrix(predO$predictions)
  dimnames(survO) <- NULL
  
  N <- dim(survO)[1]
  P <- dim(survO)[2]
  
  hP <- matrix(-99, ncol=length(times), nrow=N)
  
  for (i in 1:N) # @thomas : merci de voir si tu augmenter la vitesse du calcul de hP
  {
    hP[i,] <- sapply(times, FUN="exphaz", age=data[i,age],
                     sex=data[i,sex], year=data[i,year]) * inter
  }
  
  hcumO <- -1*log(survO)
  hinstO <- hcumO[,2:length(times)] - hcumO[,1:(length(times)-1)]
  hinstO[hinstO==Inf] <- NA
  
  for (i in 1:N)
  {
    if(sum(survO[i,]==0)>0)
    {
      hinstO[i,is.na(hinstO[i,])] <- hinstO[i,!is.na(hinstO[i,])][sum(!is.na(hinstO[i,]))]
    }
  }
  
  distOa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstO )))))
  distOb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstO )))))
  distO <- cbind(rep(0, N), (distOa + distOb)/2)
  
  hinstP <- hP[,1:(length(times)-1)]
  distPa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstP )))))
  distPb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstP )))))
  distP <- cbind(rep(0, N), (distPa + distPb)/2)
  
  hinstE <- hinstO - hinstP
  distEa <- t(as.matrix(cumsum(data.frame(t(survO[,-P] * hinstE )))))
  distEb <- t(as.matrix(cumsum(data.frame(t(survO[,-1] * hinstE )))))
  distE <- cbind(rep(0, N), (distEa + distEb)/2)
  
  distP <- distP * (1-survO)/distO
  distE <- distE * (1-survO)/distO
  
  distP[survO==1] <- 0
  distE[survO==1] <- 0
  
  distO <- distP + distE
  
  distPinf <- distP[,P]
  distEinf <- distE[,P]
  
  estimPcure <- (round(distPinf + distEinf, 10) == 1)
  
  survP <- cbind(rep(1, N), exp(-t(as.matrix(cumsum(data.frame(t(hinstP)))))))
  survU <- cbind(rep(1, N), exp(-t(as.matrix(cumsum(data.frame(t(hinstE)))))))
  
  Pcure <- distPinf / (distPinf + (1-distPinf) * survU)
  
  res <- list(formula = formula,
              data = data,
              ratetable = ratetable,
              age = age,
              sex= sex,
              year = year,
              pro.time = pro.time,
              inter = splann$inter,
              size = splann$size,
              decay = splann$decay,
              fitsurvivalnet = splann,
              predictions = list(times=times, Fc=distE, Fp = distP, tPcure = Pcure, aPcure = distPinf)
  )
  
  class(res) <- "rsPLANN"
  return(res)
}

#data(dataK)
#data(fr.ratetable)


#sp <- sPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker, data = dataK,
#                             pro.time = 365.241*20, inter=365.241/12, size = 32, decay = 0.01,
#                            maxit = 500, MaxNWts=10000)

#sp$intervals

#predO <- predict(sp, newtimes=sp$intervals)

#rs1 <- rsPLANN(Surv(time, event) ~ stade + delay + sex +  biomarker, data = dataK,
#               pro.time = 365.241*20, inter=365.241/12, size = 32, decay = 0.01,
#               maxit = 500, MaxNWts=10000, ratetable=fr.ratetable,
#               age="age", sex="sexchara", year="year")

#mdistE <- apply(rs1$predictions$Fc, FUN=mean, MARGIN=2)
#mdistP <- apply(rs1$predictions$Fp, FUN=mean, MARGIN=2)

#CPD <- cmp.rel(Surv(time, event) ~ 1,
#               data = dataK, ratetable = fr.ratetable,
#               rmap=list(age = age, sex = sex, year = year))

#plot(CPD, xlab="Post-diagnostic time (years)", xscale = 365.241, xlim=c(0,13),
#     ylab="Cumulative incidence function", ylim=c(0,1),
#     conf.int=c(1,2), col="red")
#lines(rs1$predictions$times/365.241, mdistE, type="s", col="blue")
#lines(rs1$predictions$times/365.241, mdistP, type="s", col="blue")


#plot(survfit(Surv(time, event) ~ 1, data = dataK), 
#     ylab="Patient survival", xlab="Post-diagnostic time (years)")
#lines(rs1$predictions$times, 1-(mdistE+mdistP), type="s", col="blue")


#mPcure <- apply(rs1$predictions$tPcure, FUN=mean, MARGIN=2)

#plot(rs1$predictions$times/365.241, mPcure, type="s", col="blue", 
#     xlab="Post-diagnostic time (years)", ylab="Probability of being cured")
#abline(h=0.95, col="gray", lty=2)
#abline(v=min(rs1$predictions$times[mPcure>0.95])/365.241, col="gray", lty=2)



