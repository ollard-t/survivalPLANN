
survivalPLANN <- function(formula, data, inter, size = 32, decay = 0.01,
                          maxit =100, MaxNWts=10000, ...)
{
  
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(inter)) stop("an inter argument is required")
  if (class(formula) != "formula") stop("The first argument must be a formula")
  if (class(data) != "data.frame") stop("The second argument must be a data frame")
  if (class(inter) != "numeric") stop("The inter argument must be numeric")
  
  ####### data management
  
  timevar <- as.character(formula[[2]][2])
  time <- data[,timevar]
  event <- data[,as.character(formula[[2]][3])]
  cova <- as.matrix(data[,gsub("\\+", "", attr(terms(formula), "term.labels"))])
  for(i in colnames(cova)){ if(!is.numeric(cova[,i])) stop("All covariates must be numeric")  }
  
  d <- cbind(time, event, cova)
  na <- !is.na(apply(d, MARGIN=1, FUN = "sum"))
  
  time <- time[na]
  event <- event[na]
  if (dim(cova)[2] == 1) {
    cova <- as.matrix(cova[na, ])  
  } else {
    cova <- cova[na, ]  
  }
  
  #### making of intervals from PLANN
  
  intervals <- unique(c(0, seq(inter,max(time),by=inter), max(time)))
  
  data <- data[na,]
  
  #in which interval the time t is :  
  data$Intervals = findInterval(time, intervals, left.open = TRUE)
  
  data_dup <- data[rep(seq_len(nrow(data)), data$Intervals), ]
  data_dup$Intervals <- ave(1:nrow(data_dup), data_dup[as.character(formula[[2]][2])], FUN = seq_along)
  
  
  modify_column = function(column) {
    if (length(column)==1){
      column = column
    }
    else{
      if (any(column == 1)) {  
        last_one_index = max(which(column == 1))  
        column[1:(last_one_index-1)] = 0  
      }
    }    
    return(column)
  }
  
  data_dup[,as.character(formula[[2]][3])] = ave(data_dup[,as.character(formula[[2]][3])],
                                                 floor(as.numeric(rownames(data_dup))), 
                                                 FUN = modify_column)

   #### neural network
  
  if( length( attr( terms(formula), "term.labels" ) )  == 0){
    formulaInt = as.formula(paste(as.character(formula[[2]][3]),"~",
                                  " Intervals"))
  }
else{
  formulaInt = as.formula(paste(as.character(formula[[2]][3]),"~", 
                          paste(attr(terms(formula), "term.labels"),
                          collapse = " + "),
                          "+ Intervals"))
}
  args <- list(...)
  if ("weights" %in% names(args)) {
    weights = args$weights
    dup_count <- as.numeric(table( floor(as.numeric(rownames(data_dup))) ) )
    weights <- rep(weights,dup_count)
    args$weights = weights
    print(length(args$weights))
  }

  survnet <- do.call(nnet, c(list(formula = formulaInt, data = data_dup, size = size, maxit = maxit, 
                                  MaxNWts = MaxNWts, decay = decay, entropy = TRUE), args))
  
  res <- list(formula = formula,
              fitsurvivalnet = survnet,
              data = data,
              data_dup = data_dup,
              call = survnet$call,
              inter = inter,
              size = size,
              decay = decay,
              maxit = maxit,
              MaxNWts = MaxNWts,
              coefnames = gsub("\\+", "", attr(terms(formula), "term.labels")),
              y = cbind(time = time, status = event),
              x = cova,
              intervals = intervals,
              missing = !na
              )
  class(res) <- "survivalPLANN"
  return(res)
}