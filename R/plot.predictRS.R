
plot.predictRS <- function(x, n.groups=5, type = "relative", pro.time=NULL, newdata=NULL, ...){
  
  if(!(type %in% c("relative","net","CIF_C", "CIF_P")))  stop("Argument 
                  'type' must be 'relative', 'net', 'CIF_C' or 'CIF_P' ")
  
  method_name <- NULL
  
  if(is.null(pro.time)){ pro.time <- median(x$y[,1]) }
  
  if(is.null(newdata)){
    
    cova <-data.frame(x$x)
    time <- x$y[,1];  event <- x$y[,2]
    .time <- names(x$y[1])
    .event <- names(x$y[2])
    ratetable <- x$ratetable
    
    covnames <- attr(terms(x$nnet$formula), "term.labels")
    
    ays <- x$ays
    .age <- names(x$ays)[1]
    .year <- names(x$ays)[2]
    .sex <- names(x$ays)[3]
    
    newdata <- cbind(time, event, cova, ays)
    colnames(newdata) <- c(.time, .event, covnames, .age, .year, .sex)
  }else{
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    
    ratetable <- x$ratetable
    covnames <- attr(terms(x$nnet$formula), "term.labels")
    
    .time <- names(x$y[1])
    .event <- names(x$y[2])
    .age <- names(x$ays)[1]
    .year <- names(x$ays)[2]
    .sex <- names(x$ays)[3]

    indic <- c(as.character(x$formula[[2]][2]), as.character(x$formula[[2]][3]),
               covnames,.age, .year, .sex) %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    
    cova <- data.frame(newdata[,covnames])
    time <- newdata[,as.character(x$nnet$formula[[2]][2])]
    event <- newdata[,as.character(x$nnet$formula[[2]][3])]
    
    ays = newdata[, c(.age, .year, .sex)]
    
    newdata <- cbind(time, event, cova, ays)
    colnames(newdata) <- c(.time, .event, covnames, .age, .year, .sex)
  }
 
  if( type == "relative"){
  
    .pred <- predictRS(x$nnet, data = newdata, newtimes = pro.time, 
                       ratetable = ratetable, age = .age, year = .year, sex = .sex)$ipredictions$relative_survival[,2]
    
    
  method_name = "ederer1"
  
  yname <- "rsPLANN predictions"
  xname <- "Ederer estimator predictions"
  } ##fin relative
  
  if(type == "net"){
    
    .pred <- predict(x$nnet, newdata=newdata, newtimes=pro.time,
                     ratetable = ratetable, age = .age, year = .year, sex = .sex)$ipredictions
    
    .pred <- .pred$survival_E[,1] ### on n'a plus de survie indiviudelle nette 
    #(la suvie nette étant définie comme le ratio des moyennes) 
    
    method_name = "pohar-perme"
    
    yname <- "rsPLANN predictions"
    xname <- "Pohar-Perme estimator predictions"
  }
  
  if(type == "CIF_C"){
    
    .pred <- predict(x$nnet, newdata=newdata, newtimes=pro.time,
                     ratetable = ratetable, age = .age, year = .year, sex = .sex)$ipredictions
    
    .pred <- .pred$excess_cif[,2]
    
    yname <- "rsPLANN predictions"
    xname <- "Cause specific CIF predictions"
  }
  
  if(type == "CIF_P"){
    
    .pred <- predictRS(x$nnet, data=newdata, newtimes=pro.time,
                     ratetable = ratetable, age = .age, year = .year, sex = .sex)$ipredictions
    
    .pred <- .pred$population_cif[,2]
    
    yname <- "rsPLANN predictions"
    xname <- "Population CIF predictions"
  }
  
  .grps <- as.numeric(cut(.pred,
                          breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                          labels = 1:n.groups))
  
  .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
  age = ays$age
  year =ays$year
  sex = ays$sex
  .data <- data.frame(time = time, event = event, grps = .grps, age = age, year = year, sex = sex )

  if(!is.null(method_name)){
      .obs <- c()
      .lower <- c()
      .upper <- c()
      for(i in 1:n.groups){
        datahold <- .data[.data$grps == i, ]
        
        .survfit <- summary(rs.surv(Surv(time, event) ~ 1, data = datahold,
                                ratetable = ratetable, method = method_name,
                                rmap = list(age = age, sex = sex, year = year),
                                add.times = pro.time), times = pro.time)
        
        .obs <- c(.obs, .survfit$surv)
        .lower <- c(.lower, .survfit$lower)
        .upper <- c(.upper, .survfit$upper)
        
      }
  
   }
  
  if(is.null(method_name)){
    
    
    times_list <- lapply(sort(unique(.data$grps)), function(i) .data$time[.data$event == 1 & 
                                                                            .data$grps == i ] ) 
    
    names(times_list) <- paste(".times", seq_along(times_list), sep = "")
    
    
    .all_survfit <- list()
    .all_var <- list()
    .time_fit <- c()
    .strata <- c()
    .surv <- c()
    .var <- c()
    
    
    for(i in 1:length(times_list) ){
    .survfit <- summary(cmp.rel(Surv(time, event) ~ grps, data = .data,
                                ratetable = ratetable, rmap = list(
                                  age = age, sex = sex, year = year)), times = times_list[[i]] ,
                                  scale = 1, area = FALSE)
    
    .all_survfit[[i]] <- .survfit$est[(i*2-1):(i*2),]
    .all_var[[i]] <- .survfit$var[(i*2-1):(i*2),] 
    .time_fit <- c(.time_fit, sort(times_list[[i]]))
    .strata <- c(.strata, rep(i, dim(as.data.frame(.all_survfit[[i]]))[2] ))

    
    }

    for (j in seq_along(.all_survfit)) {
      matrix_data <- as.matrix(unname(as.data.frame(.all_survfit[[j]])))
      rownames(matrix_data) <- NULL
      matrix_var <- as.matrix(unname(as.data.frame(.all_var[[j]]))) 
      rownames(matrix_var) <- NULL
      
      if (type == "CIF_C") {
        .surv <- c(.surv, matrix_data[1, ])
        .var <- c(.var, matrix_var[1, ])
      } else if (type == "CIF_P") {
        .surv <- c(.surv, matrix_data[2, ])
        .var <- c(.var,  matrix_var[2, ])
      }
    }
    
    
    .lower <- .surv - 1.96*sqrt(.var)
    .upper <- .surv + 1.96*sqrt(.var)
        
    .survfit <- list()
    .survfit$time <- .time_fit
    .survfit$strata <- .strata
    .survfit$surv <- .surv
    .survfit$lower <- .lower
    .survfit$upper <- .upper
    
  }

  # .obs <- sapply(1:n.groups, FUN = function(x) {
  #   .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
  #   .survfit$surv[ .indic ] } )
  # 
  # .lower <- sapply(1:n.groups, FUN = function(x) {
  #   .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
  #   .survfit$lower[ .indic ] } )
  # 
  # .upper <- sapply(1:n.groups, FUN = function(x) {
  #   .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
  #   .survfit$upper[ .indic ] } )
  
  
  if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
  if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
  if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
  if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
  if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
  if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
  if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
  if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
  if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
  
  if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
  if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
  
  if(hasArg(ylab)==FALSE) {ylab <- yname} else {ylab <- list(...)$ylab}
  if(hasArg(xlab)==FALSE) {xlab <- xname} else {xlab <- list(...)$xlab}
  if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
  
  plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
       type = type, col = col, lty = lty, lwd = lwd, main=main,
       pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
  
  abline(c(0,1), lty=2)
  
  segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
}