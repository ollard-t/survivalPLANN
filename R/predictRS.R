
predictRS <- function(object, data, newtimes = NULL, ratetable, age, year, sex)
{

  if (missing(object)) stop("an object of the class sPLANN is required")
  if (missing(data)) stop("a data argument is required")

  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
      splann <- object
  
      predO <- predict(object=splann, newdata=data, newtimes=splann$intervals)
      
      times <- predO$times
      
      survO <- as.matrix(predO$predictions)
      dimnames(survO) <- NULL
      
      N <- dim(survO)[1]
      P <- dim(survO)[2]
      
      hcumO <- -1*log(survO)
      hinstO <- (hcumO[,2:length(times)] - hcumO[,1:(length(times)-1)])/splann[["inter"]]
      hinstO[hinstO==Inf] <- NA
      
      for (i in 1:N)
      { 
        if(sum(survO[i,]==0)>0)
        {
          hinstO[i,is.na(hinstO[i,])] <- hinstO[i,!is.na(hinstO[i,])][sum(!is.na(hinstO[i,]))]
        }
      }
      
      hP <- matrix(-99, ncol=length(times), nrow=N)
      
      exphaz <- function(x,age,sex,year) {
        survivalNET::expectedhaz(ratetable, age=age, sex=sex, year=year, time=x)
        }
      
      for (i in 1:N) # @Thomas : merci de voir si tu augmenter la vitesse du calcul de hP
      {
        hP[i,] <- sapply(times, FUN="exphaz", age=data[i,age], sex=data[i,sex], year=data[i,year])
      }
      
      hinstP <- hP[,1:(length(times)-1)]
      
      hinstE <- hinstO - hinstP
      
      .temp <- rep(-99, N * length(1:max(times)))
      
      results <- data.frame(id=sort(rep(1:N, length(1:max(times)))), times = .temp,
                            overall_hazard=.temp, population_hazard=.temp, relative_hazard=.temp,
                            overall_survival=.temp, population_survival=.temp, relative_survival=.temp,
                            population_cif=.temp, excess_cif=.temp)
      
      for (i in 1:N)
      {
        temp0 <- data.frame(hinstO = hinstO[i, ], hinstE = hinstE[i, ], hinstP = hinstP[i, ],
                            times0 = times[-P], times1 = times[-1], interval=1:length(times[-P]))
        
        find.it <- function(x) { return(temp0$interval[x>=temp0$times0 & x<temp0$times1]) }
        
        temp1 <- data.frame(times = 1:max(times))
        
        temp1$interval <- sapply(temp1$times, FUN="find.it")
        
        temp2 <- merge(temp1, temp0, by="interval")
        
        temp2$overall_survival<- exp(-cumsum(temp2$hinstO))
        temp2$population_survival<- exp(-cumsum(temp2$hinstP))
        temp2$relative_survival <- exp(-cumsum(temp2$hinstE))
        
        temp2$overall_hazard <- temp2$hinstO
        temp2$population_hazard <- temp2$hinstP
        temp2$relative_hazard <- temp2$hinstE
        
        temp2$population_cif <- cumsum(temp2$overall_survival * temp2$hinstP) # p464 - subsection 2.1 - 2nd equation (Mozumder et al. 2017)
        temp2$excess_cif <- cumsum(temp2$overall_survival * temp2$hinstE) # p464 - subsection 2.1 - 2nd equation (Mozumder et al. 2017)
        
        results[results$id==i, -1] <- temp2[, c("times",
                                                "overall_hazard", "population_hazard", "relative_hazard",
                                                "overall_survival", "population_survival", "relative_survival",
                                                "population_cif", "excess_cif")]
      }
      
      ipredictions <- list(
        overall_survival = matrix(results$overall_survival, ncol=length(1:max(times)), byrow = TRUE),
        overall_hazard =  matrix(results$overall_hazard, ncol=length(1:max(times)), byrow = TRUE),
        population_survival = matrix(results$population_survival, ncol=length(1:max(times)), byrow = TRUE),
        population_hazard = matrix(results$population_hazard, ncol=length(1:max(times)), byrow = TRUE),
        relative_survival = matrix(results$relative_survival, ncol=length(1:max(times)), byrow = TRUE),
        relative_hazard = matrix(results$relative_hazard, ncol=length(1:max(times)), byrow = TRUE),
        population_cif = matrix(results$population_cif, ncol=length(1:max(times)), byrow = TRUE), 
        excess_cif = matrix(results$excess_cif, ncol=length(1:max(times)), byrow = TRUE)
      )
      
      .numerator <- apply(ipredictions$overall_survival, FUN="sum", MARGIN=2)
      
      # equation 6 in Biometrics (2012)
      observable_net_hazard <- apply(ipredictions$overall_survival * ipredictions$relative_hazard, FUN="sum", MARGIN=2) / .numerator
      observable_net_survival <- exp(-cumsum(observable_net_hazard))
      
      # equation 6 in Biometrics (2012) & lambda_P in JSS (2018)
      population_hazard  <- apply(ipredictions$overall_survival * ipredictions$population_hazard, FUN="sum", MARGIN=2) / .numerator
      population_survival <- exp(-cumsum(population_hazard))
      
      # from the two previous equations (just after the equation 6 in Biometrics (2012))
      overall_hazard <- apply(ipredictions$overall_survival * ipredictions$overall_hazard, FUN="sum", MARGIN=2) / .numerator
      overall_survival <- exp(-cumsum(overall_hazard))
      
      .numerator <- apply(ipredictions$relative_survival, FUN="sum", MARGIN=2)
      net_hazard <- apply(ipredictions$relative_survival * ipredictions$relative_hazard, FUN="sum", MARGIN=2) / .numerator
      net_survival <- exp(-cumsum(net_hazard))
      
      .numerator <- apply(ipredictions$population_survival, FUN="sum", MARGIN=2)
      relative_ratio_hazard <- overall_hazard - 
        apply(ipredictions$population_survival * ipredictions$population_hazard, FUN="sum", MARGIN=2)  /  .numerator
      relative_ratio_survival <- exp(-cumsum(relative_ratio_hazard))
      
      excess_cif <- cumsum( overall_survival * observable_net_hazard) # equation end page 4 in JSS (2018)
      population_cif <- cumsum( overall_survival * population_hazard) # equation end page 4 in JSS (2018)
      
      #population_cif[survO==1] <- 0
      #excess_cif[survO==1] <- 0
      
      last_population_cif <- population_cif[length(population_cif)]
      last_excess_cif <- excess_cif[length(excess_cif)]
      
      #estimPcure <- (round(distPinf + distEinf, 3) == 1)
      #Pcure <- distPinf / survU
      #Pcure[estimPcure==FALSE,] <- NA
      
      #sumSp <- apply(survP, FUN="sum", MARGIN=2)
      
      #SPlP <- survP[,1:(P-1)] * hinstP
      #sumSPlP <- apply(SlP, FUN="sum", MARGIN=2)
      
      
      if(!is.null(newtimes)) 
        {
        newtimes <- unique(newtimes)      
        
        if(0 %in% newtimes){
          newtimes <- sort(newtimes[-(newtimes == 0)])
          warning("To assure stability in the function, 0 was removed from the newtimes.")
        }
        idx <- findInterval(newtimes, times, left.open = TRUE)
        
        ipredictions$overall_survival <- as.data.frame( ipredictions$overall_survival[,pmin(idx,length(times-1))] )
        ipredictions$overall_hazard <- as.data.frame( ipredictions$overall_hazard[,pmin(idx,length(times-1))] )
        ipredictions$population_survival <- as.data.frame( ipredictions$population_survival[,pmin(idx,length(times-1))] )
        ipredictions$population_hazard <- as.data.frame( ipredictions$population_hazard[,pmin(idx,length(times-1))] )
        ipredictions$relative_survival <- as.data.frame( ipredictions$relative_survival[,pmin(idx,length(times-1))] )
        ipredictions$relative_hazard <- as.data.frame( ipredictions$relative_hazard[,pmin(idx,length(times-1))] )
        ipredictions$population_cif <- as.data.frame( ipredictions$population_cif[,pmin(idx,length(times-1))] )
        ipredictions$excess_cif <- as.data.frame( ipredictions$excess_cif[,pmin(idx,length(times-1))] )
        
        overall_survival <- overall_survival[pmin(idx,length(times-1))] 
        overall_hazard <- overall_hazard[pmin(idx,length(times-1))] 
        population_survival <- population_survival[pmin(idx,length(times-1))] 
        population_hazard <- population_hazard[pmin(idx,length(times-1))] 
        observable_net_hazard <- observable_net_hazard[pmin(idx,length(times-1))] 
        observable_net_survival <- observable_net_survival[pmin(idx,length(times-1))] 
        relative_ratio_hazard <- relative_ratio_hazard[pmin(idx,length(times-1))] 
        relative_ratio_survival <- relative_ratio_survival[pmin(idx,length(times-1))] 
        net_hazard <- net_hazard[pmin(idx,length(times-1))] 
        net_survival <- net_survival[pmin(idx,length(times-1))] 
        excess_cif <- excess_cif[pmin(idx,length(times-1))] 
        population_cif <- population_cif[pmin(idx,length(times-1))] 
       
        
        times <- newtimes
        
      }
      
      
      
      res <- list(
        nnet = splann,
        times = 0:max(times),
        ays = splann$data[,c(age, year, sex)],
        ratetable = ratetable,
        #    max_cif = list(asymptotic = estimPcure,
        #                   population = distPinf,
        #                   excess = distEinf),
        ipredictions = list(
          overall_survival = cbind(rep(1, N), ipredictions$overall_survival),
          overall_hazard =  cbind(ipredictions$overall_hazard, rep(NA, N)),
          population_survival = cbind(rep(1, N), ipredictions$population_survival),
          population_hazard = cbind(ipredictions$population_hazard, rep(NA, N)),
          relative_survival = cbind(rep(1, N), ipredictions$relative_survival),
          relative_hazard = cbind(ipredictions$relative_hazard, rep(NA, N)),
          population_cif = cbind(rep(0, N), ipredictions$population_cif), 
          excess_cif = cbind(rep(0, N), ipredictions$excess_cif)
        ),
        mpredictions = list(
          overall_survival = c(1, overall_survival),
          overall_hazard = c(overall_hazard, NA),
          population_survival = c(1, population_survival),  
          population_hazard = c(population_hazard, NA),  # equation 6 in Biometrics (2012) & lambda_P in JSS (2018)
          observable_net_hazard = c(observable_net_hazard, NA),  # equation 6 in (Biometrics (2012)
          observable_net_survival = c(1, observable_net_survival),
          relative_ratio_hazard = c(relative_ratio_hazard, NA), # equation 9 in Biometrics (2012)
          relative_ratio_survival = c(1, relative_ratio_survival),
          net_hazard = c(net_hazard, NA), # equation 4 in Biometrics (2012)
          net_survival = c(1, net_survival), # net survival from equation 4 in Biometrics (2012)
          excess_cif =  c(0, excess_cif), # equation page 4 in JSS (2018)
          population_cif =  c(0, population_cif) # equation page 4 in JSS (2018)
          #      cure = apply(Pcure, FUN="mean", MARGIN=2)
        )
      )
      class(res) <- "predictRS"
      return(res)
}