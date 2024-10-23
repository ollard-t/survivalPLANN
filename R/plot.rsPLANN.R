
plot.rsPLANN <- function(x, n.groups=5, pro.time=NULL, newdata=NULL,
                         ratetable, ...){
  
  if(is.null(pro.time)){ pro.time <- median(x$y[,1]) }
  
  if(is.null(newdata))
  {
    cova <-data.frame(x$x)
    time <- x$y[,1];  event <- x$y[,2]
    .age <- x$ays$age; .year <- x$ays$year; .sex <- x$ays$sex
  }
  
  if(!is.null(newdata))
  {
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    
    indic <- gsub("\\+", "", attr(terms(x$formula), "term.labels") ) %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    
    cova <- data.frame(newdata[,gsub("\\+", "", attr(terms(x$formula), "term.labels"))])
    time <- newdata[,as.character(x$formula[[2]][2])]
    event <- newdata[,as.character(x$formula[[2]][3])]
    .age <- newdata[,.age]; .sex <- newdata[,.sex]; .year <- newdata[,.year]
  }
  
}