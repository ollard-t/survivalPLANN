
predict.rsPLANN <- function(object, newdata = NULL, newtimes = NULL, ratetable, age, year, sex, ...) # @Thomas: inclure ratetable et ses arguments dans formula
{
     
#  object <- rs1
#  newdata = dataK[1:200,]
#  newtimes = NULL
#  ratetable <- fr.ratetable
#  age="age"
 # sex="sexchara"
#  year="year"
  
    if( is.null(newdata) & is.null(newtimes))
    {     res <- list(times = object$predictions$times, predictions = object$predictions[-1]) 
    } else 
      
    {
      splann <- object$fitsurvivalnet
      
      if(is.null(newdata)) {data <- splann$data} else {data <- newdata}
      
      predO <- predict(object=splann, newdata=data, newtimes=splann$intervals)
      
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
                         sex=data[i,sex], year=data[i,year]) * splann$inter
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
      
      if(!is.null(newtimes)) # @Thomas : merci de vérifier que je renvois les bonnes valeurs dans cette condition -> voir plot dans l'exemple du fichier Rd
        {
      idx <- findInterval(newtimes, times, left.open = TRUE)
      distE <- distE[,pmin(idx,length(times-1))]
      distP <- distP[,pmin(idx,length(times-1))]
      Pcure <- Pcure[,pmin(idx,length(times-1))]
      times <- newtimes
        }
      
      res <- list(times=times,
                  predictions = list(times=times, Fc=distE, Fp = distP, tPcure = Pcure, aPcure = distPinf) )
    }
  
  return(res)
}