
predict.survivalPLANN <- function(object, newdata = NULL, newtimes = NULL, ...)
{
  
  intervals = object$intervals
  formula = object$formula
  time = object$y[,1]
  
  if(!is.null(newdata))
  {
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    indic <- gsub("\\+", "", attr(terms(formula), "term.labels") ) %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    covariates <- as.matrix(newdata[, gsub("\\+", "", attr(terms(object$formula), "term.labels"))])
  }
  
  if(is.null(newdata))  { newdata <- object$data[,-dim(object$data)[2]]}


  data_dupli = newdata[rep(seq_len(nrow(newdata)), each = length(intervals)-1), ]
  data_dupli$Intervals <- ave(1:nrow(data_dupli),
                              floor(as.numeric(rownames(data_dupli))),
                              FUN = seq_along)
  
  
  predictions <- predict(object$fitsurvivalnet, data_dupli, type ="raw", ...)

  grouped_df <- split(cbind(data_dupli,predictions),
                     rep(1:ceiling(nrow(data_dupli)/(length(intervals)-1)),
                     each=length(intervals)-1, length.out=nrow(data_dupli)))


  result_group <- lapply(grouped_df, function(group) {
    group$Survie <- cumprod(1-group$predictions)
    return(group)
  })
  result_df <- do.call(rbind, result_group)$Survie

  predictions <- as.data.frame(matrix(result_df, ncol = length(intervals)-1, byrow = TRUE))
  predictions <- cbind(rep(1,dim(predictions)[1]), predictions)

  ints_names <- c()

  for( i in 1:(length(object$intervals)-1) ){
    ints_names = c(ints_names, paste0("(",round(object$intervals[i], digits = 2), ";",
                                      round(object$intervals[i+1], digits = 2),"]"))
  }
  colnames(predictions) = c("0",ints_names)
  
  if(is.null(newtimes))  { newtimes <- sort(unique(time))}
  else{ 
    if(!is.vector(newtimes))stop("newtimes must be a vector")
    if(any(max(time)<newtimes))warning("The values of 'newtimes' are  greater than the max(time) of your training data base")    
    idx=findInterval(newtimes, intervals, left.open = TRUE)
    predictions = predictions[,pmin(idx+1,length(intervals-1))]
  }
  

  res <- list(times = newtimes, predictions = predictions)
  
  return(res)
  
}