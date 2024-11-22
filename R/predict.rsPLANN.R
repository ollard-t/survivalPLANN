
predict.rsPLANN <- function(object, newdata = NULL, newtimes = NULL, ...)
{
  
  if( is.null(newdata) & is.null(newtimes))
  {     res <- list(times = object$times, ipredictions = object$ipredictions,
                    mpredictions = object$mpredictions) 
  } else 
    
  {
    splann <- object$fitsurvivalnet
    ratetable <- object$ratetable
    
    all_terms <- attr(terms(object$formula), "term.labels")
    
    covnames <- names(as.data.frame(object$fitsurvivalnet$x))
    
    ratetable_terms <- grep("ratetable\\(", all_terms, value = TRUE)
    extract_vars <- function(term) {
      var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
      vars <- trimws(unlist(strsplit(var_string, ",")))
      return(vars)
    }
    
    assign_ratetable_vars <- function(vars) {
      age <- year <- sex <- NULL
      for (var in vars) {
        if (grepl("age = ", var)) {
          age <- sub("age = ", "", var)
        } else if (grepl("year = ", var)) {
          year <- sub("year = ", "", var)
        } else if (grepl("sex = ", var)) {
          sex <- sub("sex = ", "", var)
        }
      }
      unnamed_vars <- setdiff(vars, c(age, sex, year))
      if (length(unnamed_vars) > 0) {
        if (is.null(age) && length(unnamed_vars) >= 1) age <- unnamed_vars[1]
        if (is.null(year) && length(unnamed_vars) >= 2) year <- unnamed_vars[2]
        if (is.null(sex) && length(unnamed_vars) >= 3) sex <- unnamed_vars[3]
      }
      return(list(age = age, sex = sex, year = year))
    }
    
    extract_vars <- function(term) {
      var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
      vars <- trimws(unlist(strsplit(var_string, ",")))
      return(vars)
    }
    
    ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
    
    age <- ratetable_vars$age
    year <- ratetable_vars$year
    sex <- ratetable_vars$sex
    
    if(is.null(newdata)) {
      
      data <- splann$data
    
      }else{
      
        if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
      
        indic <- c(age,year,sex, covnames)  %in% names(newdata) 
        if( sum(!indic) > 0 ) stop("Missing predictor in the data frame. 'newdata' also needs
                                   'age', 'sex' and 'year' for the ratetable.")
      
        data <- newdata
      }
    
    predO <- predict(object=splann, newdata=data, newtimes=splann$intervals)
    
    times <- predO$times
    
    exphaz <- function(x,age,sex,year) { survivalNET::expectedhaz(ratetable, age=age, sex=sex, year=year, time=x)}
    
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
    
    sumS <- apply((1-distO), FUN="sum", MARGIN=2)
    
    SlE <- (1-distO) * cbind(rep(0, N), hinstE)
    sumSlE <- apply(SlE, FUN="sum", MARGIN=2)
    lambda_C <- sumSlE/sumS
    #Lambda_C <- cumsum(lambda_C)
    
    SlP <- (1-distO) * cbind(rep(0, N), hinstP)
    sumSlP <- apply(SlP, FUN="sum", MARGIN=2)
    lambda_P <- sumSlP/sumS
    #Lambda_P <- cumsum(lambda_P)
    
    # warning -> NA pour tCure ...
    
    if(!is.null(newtimes)) # @Thomas : merci de vérifier que je renvois les bonnes valeurs dans cette condition -> voir plot dans l'exemple du fichier Rd
    {
      idx <- findInterval(newtimes, times, left.open = TRUE)
      
      distO <- as.data.frame( distO[,pmin(idx,length(times-1))] )
      distE <- as.data.frame( distE[,pmin(idx,length(times-1))] )
      distP <- as.data.frame( distP[,pmin(idx,length(times-1))] )
      Pcure <- as.data.frame( Pcure[,pmin(idx,length(times-1))] )
      survP <- as.data.frame( survP[,pmin(idx,length(times-1))] )
      survU <- as.data.frame( survU[,pmin(idx,length(times-1))] )
      times <- newtimes
      
      names(distO) <- times
      names(distE) <- times
      names(distP) <- times
      names(Pcure) <- times
      names(survP) <- times
      names(survU) <- times
      
          }
    
    res <- list(
      times = times,
      ipredictions = list(survival_O=1-distO,
                          survival_P=survP,
                          survival_R=(1-distO)/survP,
                          survival_E=survU, # remarque : S(1-distO)/survP = survU
                          CIF_C = distE, CIF_P = distP, maxCIF_P = distPinf,
                          cure = Pcure),
      mpredictions = list(survival_O = apply((1-distO), FUN="mean", MARGIN=2),
                          survival_P = apply(survP, FUN="mean", MARGIN=2),
                          survival_R = apply((1-distO), FUN="mean", MARGIN=2)/apply(survP, FUN="mean", MARGIN=2),
                          survival_E = apply((1-distO)/survP, FUN="mean", MARGIN=2),
                          CIF_C =  apply(distE, FUN="mean", MARGIN=2),
                          CIF_P =  apply(distP, FUN="mean", MARGIN=2),   
                          cure = apply(Pcure, FUN="mean", MARGIN=2),
                          hazard_P = lambda_P,
                          hazard_C = lambda_C )
    )
    
  }
  
  return(res)
}