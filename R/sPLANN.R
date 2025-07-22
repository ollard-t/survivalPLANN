
sPLANN <- function(formula, data, pro.time=NULL, inter, size = 32, decay = 0.01,
                          maxit =100, MaxNWts=10000, trace = FALSE, ...)
{
  
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(inter)) stop("an inter argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (as.character(class(inter)) != "numeric") stop("The inter argument must be numeric")
  
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
  
  data <- data[na,]
  
  #### making of intervals from PLANN
  
  # if (is.null(pro.time)) {pro.time <- max(time[event==1])} # modif
  if (is.null(pro.time)) {pro.time <- max(time)} # modif
  
  intervals <- unique(c(0, seq(inter, pro.time, by=inter), pro.time)) # modif
  
  #in which interval the time t is :  
  data$Intervals = findInterval(time, intervals, left.open = TRUE)
  
  data$id <- seq_len(nrow(data))
  data_dup <- data[rep(data$id, data$Intervals), ]
  data_dup$id <- rep(data$id, data$Intervals)
  data_dup$Intervals <- ave(data_dup$id, data_dup$id, FUN = seq_along)
  data_dup <- data_dup[,-dim(data_dup)[2]] 
 
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
  }else{ #chgt
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
  }

  survnet <- do.call(nnet, c(list(formula = formulaInt, data = data_dup, size = size, maxit = maxit, 
                                  MaxNWts = MaxNWts, decay = decay, entropy = TRUE, trace = trace), args))
  
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
              y = cbind(time, event),
              x = cova,
              intervals = intervals,
              pro.time = pro.time, #pour calcul logll predictRS
              missing = !na
              )
  class(res) <- "sPLANN"
  return(res)
}