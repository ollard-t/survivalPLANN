
predictRS <- function(object, data, newtimes = NULL, ratetable, age, year, sex)
{
  
  if (!inherits(object, "sPLANN")) stop("The object must be of class 'sPLANN'")
  if (missing(object)) stop("an object of the class sPLANN is required")
  if (missing(data)) stop("a data argument is required")

  if (missing(ratetable)) stop("a ratetable argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
  covnames <- colnames(object$x)
  .age <- age
  .year <- year
  .sex <- sex
  indic <- c(as.character(object$formula[[2]][2]), as.character(object$formula[[2]][3]),
             covnames,.age, .year, .sex) %in% names(data) 
  if( sum(!indic) > 0 ) stop("Missing predictor in the new data frame")
  
  if(!is.null(newtimes) & 0 %in% newtimes)stop("Error : 0 found in 'newtimes'. For stablity of the function, 0 can't be in 'newtimes'.")
      
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
      
      exphaz <- function(x,age,sex,year,max_age,max_year) {
        survivalNET::expectedhaz(ratetable, age=age, sex=sex, year=year, time=x, max_age = max_age, max_year= max_year)
        }
      
      max_age <- max(as.numeric(dimnames(ratetable)[[1]]))
      max_year <- max(as.numeric(dimnames(ratetable)[[2]]))
      
      for (i in 1:N) # @Thomas : merci de voir si tu augmenter la vitesse du calcul de hP
      {
        hP[i,] <- sapply(times, FUN="exphaz", age=data[i,age], sex=data[i,sex], year=data[i,year], 
                         max_age = max_age, 
                         max_year = max_year)
      }
      
      hinstP <- hP[,1:(length(times)-1)]
      
      hinstE <- pmax(hinstO - hinstP, 0)
      
      if(is.null(newtimes)){
        
        .temp <- rep(-99, N * (length(times)-1))
      
        results <- data.frame(id=sort(rep(1:N, (length(times)-1) )), times = .temp,
                            overall_hazard=.temp, population_hazard=.temp, relative_hazard=.temp,
                            overall_survival=.temp, population_survival=.temp, relative_survival=.temp,
                            population_cif=.temp, excess_cif=.temp)
      }else{
        
        .temp <- rep(-99, N * (length(times)-1+length(newtimes)))
        
        results <- data.frame(id=sort(rep(1:N, (length(times)-1+length(newtimes)) )), times = .temp,
                              overall_hazard=.temp, population_hazard=.temp, relative_hazard=.temp,
                              overall_survival=.temp, population_survival=.temp, relative_survival=.temp,
                              population_cif=.temp, excess_cif=.temp)
      }
      
      for (i in 1:N)
      { 
        temp0 <- data.frame(hinstO = hinstO[i, ], hinstE = hinstE[i, ], hinstP = hinstP[i, ],
                            times0 = times[-P], times1 = times[-1], interval=1:length(times[-P]))
        
        temp1 <- data.frame(times = times[-1]) 
        
        temp1$interval <- findInterval(temp1$times, temp0$times0, left.open = TRUE)

        # temp2 <- merge(temp1, temp0, by="interval")
        temp2 <- temp0[temp1$interval,]
        temp2$times <- temp1$times
        
        if(!is.null(newtimes)){
          row_index <- findInterval(newtimes, temp0$times0, left.open = TRUE)
          new_row <- temp2[row_index, ]
          new_row$times <- newtimes  
          temp2 <- rbind(temp2, new_row)
          temp2 <- temp2[order(temp2$times),]
        }
        mult <- diff(floor(c(0,temp2$times)))

        temp2$overall_survival<- exp(-cumsum(mult*temp2$hinstO))
        temp2$population_survival<- exp(-cumsum(mult*temp2$hinstP))
        temp2$relative_survival <- exp(-cumsum(mult*temp2$hinstE))
        
        temp2$overall_hazard <- temp2$hinstO
        temp2$population_hazard <- temp2$hinstP
        temp2$relative_hazard <- temp2$hinstE
        
        temp2$population_cif <- cumsum(mult*temp2$overall_survival * temp2$hinstP) # p464 - subsection 2.1 - 2nd equation (Mozumder et al. 2017)
        temp2$excess_cif <- cumsum(mult*temp2$overall_survival * temp2$hinstE) # p464 - subsection 2.1 - 2nd equation (Mozumder et al. 2017)
        
        results[results$id == i, -1] <- temp2[, c("times",
                                                "overall_hazard", "population_hazard", "relative_hazard",
                                                "overall_survival", "population_survival", "relative_survival",
                                                "population_cif", "excess_cif")]
      }
      
      ipredictions <- list(
        overall_survival = matrix(results$overall_survival, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        overall_hazard =  matrix(results$overall_hazard, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        population_survival = matrix(results$population_survival, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        population_hazard = matrix(results$population_hazard, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        relative_survival = matrix(results$relative_survival, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        relative_hazard = matrix(results$relative_hazard, ncol=length(times[-1])+length(newtimes), byrow = TRUE),
        population_cif = matrix(results$population_cif, ncol=length(times[-1])+length(newtimes), byrow = TRUE), 
        excess_cif = matrix(results$excess_cif, ncol=length(times[-1])+length(newtimes), byrow = TRUE)
      )
      
      .numerator <- apply(ipredictions$overall_survival, FUN="sum", MARGIN=2)
      
      # equation 6 in Biometrics (2012)
      observable_net_hazard <- apply(ipredictions$overall_survival * ipredictions$relative_hazard, FUN="sum", MARGIN=2) / .numerator
      observable_net_survival <- exp(-cumsum(mult*observable_net_hazard)) ##pred très proche de ce qu'on avait avant mais petite différence /!\
      
      # equation 6 in Biometrics (2012) & lambda_P in JSS (2018)
      population_hazard  <- apply(ipredictions$overall_survival * ipredictions$population_hazard, FUN="sum", MARGIN=2) / .numerator
      population_survival <- exp(-cumsum(mult*population_hazard)) ##pred très proche (3 chiffres après la virgule) de ce qu'on avait avant mais petite différence /!\
      
      # from the two previous equations (just after the equation 6 in Biometrics (2012))
      overall_hazard <- apply(ipredictions$overall_survival * ipredictions$overall_hazard, FUN="sum", MARGIN=2) / .numerator
      overall_survival <- exp(-cumsum(mult*overall_hazard)) ##pred très proche de ce qu'on avait avant mais petite différence /!\
      
      .numerator <- apply(ipredictions$relative_survival, FUN="sum", MARGIN=2)
      net_hazard <- apply(ipredictions$relative_survival * ipredictions$relative_hazard, FUN="sum", MARGIN=2) / .numerator # equaition 4 biometrics (2012)
      net_survival <- exp(-cumsum(mult*net_hazard)) ##pred très proche de ce qu'on avait avant mais petite différence /!\
      
      .numerator <- apply(ipredictions$population_survival, FUN="sum", MARGIN=2)
      relative_ratio_hazard <- overall_hazard - 
        apply(ipredictions$population_survival * ipredictions$population_hazard, FUN="sum", MARGIN=2)  /  .numerator
      relative_ratio_survival <- exp(-cumsum(mult*relative_ratio_hazard)) ##pred très proche de ce qu'on avait avant mais petite différence /!\
      
      excess_cif <- cumsum( mult*overall_survival * observable_net_hazard) # equation end page 4 in JSS (2018)
      ##pred très proche de ce qu'on avait avant mais petite différence /!\
      population_cif <- cumsum( mult*overall_survival * population_hazard) # equation end page 4 in JSS (2018)
      ##pred très proche de ce qu'on avait avant mais petite différence /!\
      
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
      
      ### log-likelihood 
      ## on récpère intervalles où tombent les temps d'evt et de censure
      
      event_time <- findInterval(splann$y[,1], splann$intervals,left.open = TRUE)
      
      # #on récupère le risque instantané indivudel observé au temps d'evt/censure
      # ind_hinstO <- sapply(1:(dim(data)[1]), function(i) {
      #   hinstO[i, event_time[i]]
      # })
      # #et la survie observée individuelle au temps d'evt/censure 
      # ind_survO <-  sapply(1:(dim(data)[1]), function(i) {
      #   ipredictions$overall_survival[i, event_time[i]]
      # })
      # 
      # loglik <- sum(splann$y[,2]*log(ind_hinstO)+log(ind_survO)) 
      
      ### log-likelihood 
      
      #on récupère le risque instantané populationnel au temps d'evt/censure
      pop_hinst <- sapply(1:(dim(data)[1]), function(i) {
        ipredictions$population_hazard[i, event_time[i]]
      })
      #on récupère le risque instantané en exces au temps d'evt/censure
      exc_hinst <- sapply(1:(dim(data)[1]), function(i) {
        ipredictions$relative_hazard[i, event_time[i]]
      })
      #et la survie nette au temps d'evt/censure 
      net_surv <-  sapply(1:(dim(data)[1]), function(i) {
        ipredictions$relative_survival[i, event_time[i]]
      })
      
      if (identical(data, splann$data[,-c(length(splann$data)-1 ,length(splann$data))])) {
          loglik <- sum(splann$y[,2]*log(pop_hinst+exc_hinst)+log(net_surv)) 
      }
          
      if(!is.null(newtimes)) 
        {
        newtimes <- unique(newtimes)      
        
        # if(0 %in% newtimes){
        #   newtimes <- sort(newtimes[-(newtimes == 0)])
        #   warning("To assure stability in the function, 0 was removed from the newtimes.")
        # }
        nouveautime <- sort(c(newtimes, times))
        idx <- findInterval(newtimes, nouveautime, left.open = TRUE) ## c'était times avant nouveautime (et pareil après si on doit rechanger)
        
        ipredictions$overall_survival <- as.data.frame( ipredictions$overall_survival[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$overall_hazard <- as.data.frame( ipredictions$overall_hazard[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$population_survival <- as.data.frame( ipredictions$population_survival[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$population_hazard <- as.data.frame( ipredictions$population_hazard[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$relative_survival <- as.data.frame( ipredictions$relative_survival[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$relative_hazard <- as.data.frame( ipredictions$relative_hazard[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$population_cif <- as.data.frame( ipredictions$population_cif[,pmin(idx,length(nouveautime)-1)] )
        ipredictions$excess_cif <- as.data.frame( ipredictions$excess_cif[,pmin(idx,length(nouveautime)-1)] )
        
        overall_survival <- overall_survival[pmin(idx,length(nouveautime)-1)] 
        overall_hazard <- overall_hazard[pmin(idx,length(nouveautime)-1)] 
        population_survival <- population_survival[pmin(idx,length(nouveautime)-1)] 
        population_hazard <- population_hazard[pmin(idx,length(nouveautime)-1)] 
        observable_net_hazard <- observable_net_hazard[pmin(idx,length(nouveautime)-1)] 
        observable_net_survival <- observable_net_survival[pmin(idx,length(nouveautime)-1)] 
        relative_ratio_hazard <- relative_ratio_hazard[pmin(idx,length(nouveautime)-1)] 
        relative_ratio_survival <- relative_ratio_survival[pmin(idx,length(nouveautime)-1)] 
        net_hazard <- net_hazard[pmin(idx,length(nouveautime)-1)] 
        net_survival <- net_survival[pmin(idx,length(nouveautime)-1)] 
        excess_cif <- excess_cif[pmin(idx,length(nouveautime)-1)] 
        population_cif <- population_cif[pmin(idx,length(nouveautime)-1)] 
       
        
        times <- newtimes
        
      }
      
      ## pour nommer les colonnes des estimations individuelles 
      ints_names <- c()
      
      for( i in 1:(length(object$intervals)-1) ){
        ints_names <- c(ints_names, paste0("(",round(object$intervals[i], digits = 2), ";",
                                           round(object$intervals[i+1], digits = 2),"]"))
      }
      
      usable_times <- times[times != 0]
      idx <- findInterval(usable_times, object$intervals, left.open = TRUE)
      col_names <- paste0(usable_times," in ",ints_names[idx])
      
      
      
      res <- list(
        nnet = splann,
        times = 0:max(times),
        x = data[,c(colnames(splann$x))],
        y =  data[, c(as.character(splann$formula[[2]][2]), as.character(splann$formula[[2]][3]))],
        ays = data[,c(age, year, sex)],
        ratetable = ratetable,
        #    max_cif = list(asymptotic = estimPcure,
        #                   population = distPinf,
        #                   excess = distEinf),
        ipredictions = list(
          overall_survival = { mat <- cbind(rep(1, N), ipredictions$overall_survival)
          colnames(mat) <- c(0,col_names)
          mat },
          overall_hazard =  { mat <- cbind(ipredictions$overall_hazard, rep(NA, N))
          colnames(mat) <- c(col_names, "NA")
          mat},
          population_survival = { mat <- cbind(rep(1, N), ipredictions$population_survival)
          colnames(mat) <- c(0,col_names)
          mat },
          population_hazard = { mat <- cbind(ipredictions$population_hazard, rep(NA, N))
          colnames(mat) <- c(col_names, "NA")
          mat},
          relative_survival = { mat <-cbind(rep(1, N), ipredictions$relative_survival)
          colnames(mat) <- c(0,col_names)
          mat },
          relative_hazard = { mat <- cbind(ipredictions$relative_hazard, rep(NA, N))
          colnames(mat) <- c(col_names, "NA")
          mat},
          population_cif = { mat <- cbind(rep(0, N), ipredictions$population_cif)
          colnames(mat) <- c(0,col_names)
          mat }, 
          excess_cif = { mat <- cbind(rep(0, N), ipredictions$excess_cif)
          colnames(mat) <- c(0,col_names)
          mat }
        ),
        mpredictions = list(
          overall_survival =  c(1, overall_survival),
          overall_hazard =  c(overall_hazard, NA),
          population_survival =  c(1, population_survival),  
          population_hazard =  c(population_hazard, NA),  # equation 6 in Biometrics (2012) & lambda_P in JSS (2018)
          observable_net_hazard =  c(observable_net_hazard, NA),  # equation 6 in (Biometrics (2012)
          observable_net_survival =  c(1, observable_net_survival),
          relative_ratio_hazard =  c(relative_ratio_hazard, NA), # equation 9 in Biometrics (2012)
          relative_ratio_survival = c(1, relative_ratio_survival),
          net_hazard =  c(net_hazard, NA), # equation 4 in Biometrics (2012)
          net_survival = c(1, net_survival), # net survival from equation 4 in Biometrics (2012)
          excess_cif =  c(0, excess_cif), # equation page 4 in JSS (2018)
          population_cif =  c(0, population_cif) # equation page 4 in JSS (2018)
          #      cure = apply(Pcure, FUN="mean", MARGIN=2)
        )
      )
      if (exists("loglik")) {
        res$loglik = loglik
      }
      class(res) <- "predictRS"
      return(res)
}