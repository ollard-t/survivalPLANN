
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
  indic <- c(covnames,.age, .year, .sex) %in% names(data) 
  if( sum(!indic) > 0 ) stop("Missing predictor in the new data frame")
  
  if(!is.null(newtimes) & 0 %in% newtimes)stop("Error : 0 found in 'newtimes'. For stablity of the function, 0 can't be in 'newtimes'.")
      
      splann <- object
      
      times <- splann$intervals
      if(!is.null(newtimes)){
        matches <- sapply(newtimes, function(nt) any(abs(times - nt) < 1e-6)) ## pour n'ajouter que des temps nouveaux dans la base
        times <- sort(c(times, newtimes[!matches]))}
      
      predO <- predict(object=splann, newdata=data, newtimes=times)
      
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
      
      max_age <- max(as.numeric(dimnames(ratetable)[[1]]))
      max_year <- max(as.numeric(dimnames(ratetable)[[2]]))
      
      hP <- t(vapply(seq_len(nrow(data)),function(i)
        survivalNET::expectedhaz(
            ratetable, age=data[i,age], sex=data[i,sex], year=data[i,year] ,
            time = times, max_age = max_age, max_year = max_year
          ),
        FUN.VALUE = numeric(length(times))
      ))
      
      
      hinstP <- hP[,1:(length(times)-1)]
      
      hinstE <- pmax(hinstO - hinstP, 0)
      
      # if(is.null(newtimes)){
        
        .temp <- rep(-99, N * (length(times)-1))
      
        results <- data.frame(id=sort(rep(1:N, (length(times)-1) )), times = .temp,
                            overall_hazard=.temp, population_hazard=.temp, relative_hazard=.temp,
                            overall_survival=.temp, population_survival=.temp, relative_survival=.temp,
                            population_cif=.temp, excess_cif=.temp)
      # }else{
      #   
      #   .temp <- rep(-99, N * (length(times)-1))
      #   
      #   results <- data.frame(id=sort(rep(1:N, (length(times)-1) )), times = .temp,
      #                         overall_hazard=.temp, population_hazard=.temp, relative_hazard=.temp,
      #                         overall_survival=.temp, population_survival=.temp, relative_survival=.temp,
      #                         population_cif=.temp, excess_cif=.temp)
      # }
      
        interval_index <- findInterval(times[-1], times[-P], left.open = TRUE)
        mult <- diff(floor(c(0, times[-1])))
        Ti <- length(times) - 1
        start_idx <- ((seq_len(N) - 1) * Ti) + 1
        end_idx   <- start_idx + Ti - 1
        
        res_mat <- matrix(NA_real_, nrow = N * Ti, ncol = 9)
        colnames(res_mat) <- colnames(results)[-1]
        
        for (i in 1:N){ 
          
          rows <- start_idx[i]:end_idx[i]
          
          hO <- hinstO[i, interval_index]
          hP <- hinstP[i, interval_index]
          hE <- hinstE[i, interval_index]
          
          cumO <- cumsum(mult * hO)
          cumP <- cumsum(mult * hP)
          cumE <- cumsum(mult * hE)
          
          overall_survival    <- exp(-cumO)
          population_survival <- exp(-cumP)
          relative_survival   <- exp(-cumE)
          
          population_cif <- cumsum(mult * overall_survival * hP)
          excess_cif     <- cumsum(mult * overall_survival * hE)
          
          res_mat[rows, ] <- cbind(
            times[-1],
            hO, hP, hE,
            overall_survival,
            population_survival,
            relative_survival,
            population_cif,
            excess_cif
          )
        }
        results[-1] <- as.data.frame(res_mat)
      
      ipredictions <- list(
        overall_survival = matrix(results$overall_survival, ncol=length(times[-1]), byrow = TRUE),
        overall_hazard =  matrix(results$overall_hazard, ncol=length(times[-1]), byrow = TRUE),
        population_survival = matrix(results$population_survival, ncol=length(times[-1]), byrow = TRUE),
        population_hazard = matrix(results$population_hazard, ncol=length(times[-1]), byrow = TRUE),
        relative_survival = matrix(results$relative_survival, ncol=length(times[-1]), byrow = TRUE),
        relative_hazard = matrix(results$relative_hazard, ncol=length(times[-1]), byrow = TRUE),
        population_cif = matrix(results$population_cif, ncol=length(times[-1]), byrow = TRUE), 
        excess_cif = matrix(results$excess_cif, ncol=length(times[-1]), byrow = TRUE)
      )
      
      .numerator <- colSums(ipredictions$overall_survival)
      
      # equation 6 in Biometrics (2012)
      observable_net_hazard <- colSums(ipredictions$overall_survival * ipredictions$relative_hazard) / .numerator
      observable_net_survival <- exp(-cumsum(mult*observable_net_hazard)) 
      
      # equation 6 in Biometrics (2012) & lambda_P in JSS (2018)
      population_hazard  <- colSums(ipredictions$overall_survival * ipredictions$population_hazard) / .numerator
      population_survival <- exp(-cumsum(mult*population_hazard)) 
      
      # from the two previous equations (just after the equation 6 in Biometrics (2012))
      overall_hazard <- colSums(ipredictions$overall_survival * ipredictions$overall_hazard) / .numerator
      overall_survival <- exp(-cumsum(mult*overall_hazard)) 
      
      .numerator <- colSums(ipredictions$relative_survival)
      net_hazard <- colSums(ipredictions$relative_survival * ipredictions$relative_hazard) / .numerator # equaition 4 biometrics (2012)
      net_survival <- exp(-cumsum(mult*net_hazard)) 
      
      .numerator <- colSums(ipredictions$population_survival)
      relative_ratio_hazard <- overall_hazard - 
        colSums(ipredictions$population_survival * ipredictions$population_hazard)  /  .numerator
      relative_ratio_survival <- exp(-cumsum(mult*relative_ratio_hazard)) 
      
      excess_cif <- cumsum( mult*overall_survival * observable_net_hazard) # equation end page 4 in JSS (2018)
      
      population_cif <- cumsum( mult*overall_survival * population_hazard) # equation end page 4 in JSS (2018)
      
      
      last_population_cif <- population_cif[length(population_cif)]
      last_excess_cif <- excess_cif[length(excess_cif)]
      
      ### log-likelihood 
      ## on récpère intervalles où tombent les temps d'evt et de censure
      if(all(c(as.character(object$formula[[2]][2]), as.character(object$formula[[2]][3])) %in% names(data))) {
        
        event_time <- findInterval(data[,as.character(object$formula[[2]][2])], splann$intervals,left.open = TRUE)
        event_time[event_time > dim(ipredictions$population_hazard)[2]] <- dim(ipredictions$population_hazard)[2]
        
        idx <- cbind(seq_len(N), event_time)
        pop_hinst <- ipredictions$population_hazard[idx]
        exc_hinst <- ipredictions$relative_hazard[idx]
        net_surv  <- ipredictions$relative_survival[idx]
        
     
        loglik <- sum(data[,as.character(object$formula[[2]][3])]*log(pop_hinst+exc_hinst)+log(net_surv))

    }
      # if(!(identical(data, splann$data[,-c(length(splann$data)-1 ,length(splann$data))]))) {
      #   if(all(c(as.character(object$formula[[2]][2]), as.character(object$formula[[2]][3])) %in% names(data))) {
      #       loglik <- sum(data[,as.character(object$formula[[2]][3])]*log(pop_hinst+exc_hinst)+log(net_surv))
      #     }else{
      #       #do nothing and don"t compute the loglikelihood
      #     }
      #   }
          
      if(!is.null(newtimes)) 
        {
        newtimes <- unique(newtimes)      
        
        nouveautime <- sort(unique(round(c(newtimes, times),6)))
        idx <- findInterval(newtimes, nouveautime, left.open = TRUE) 
      
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
        times = c(0,times), #anciennement 0:max(times) au 24 sep
        x = data[,c(colnames(splann$x))],
        ays = data[,c(age, year, sex)],
        ratetable = ratetable,
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
        )
      )
      if (exists("loglik")) {
        res$loglik = loglik
      }
      if(all(c(as.character(object$formula[[2]][2]), as.character(object$formula[[2]][3])) %in% names(data))){
        res$y =  data[, c(as.character(splann$formula[[2]][2]), as.character(splann$formula[[2]][3]))]
      }
      class(res) <- "predictRS"
      return(res)
}