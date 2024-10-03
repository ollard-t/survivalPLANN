
rsPLANN <- function(formula, data, pro.time=NULL, inter, size= 32, decay=0.01,
                    maxit=100, MaxNWts=10000, trace=FALSE,
                    ratetable) 
{
  
  ####### check errors
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(inter)) stop("an inter argument is required")
  if (as.character(class(formula)) != "formula") stop("The formula argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (as.character(class(inter)) != "numeric") stop("The inter argument must be numeric")
  
  if (missing(ratetable)) stop("a ratetable argument is required")

  
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
  ####### data management
  
  all_terms <- attr(terms(formula), "term.labels")
  ratetable_terms <- grep("ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  if(length(ratetable_terms)>1) stop("More than one 'ratetable' term found in  the formula.")
  
  formula_string <- as.character(formula)
  rhs <- gsub("\\+ ratetable\\(.*?\\)", "", formula_string[3]) 
  formula_string <- paste(formula_string[2], "~", rhs)
  formula <- as.formula(formula_string)
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
    return(list(age = age, year = year, sex = sex))
  }
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
  age <- data[,ratetable_vars$age]  
  year <- data[,ratetable_vars$year]
  sex <- data[,ratetable_vars$sex] 
  
  splann <- sPLANN(formula, data=data, pro.time=pro.time, inter=inter, 
                          size=size, decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  predO <- predict(splann, newtimes=splann$intervals)
  
  times <- predO$times
  
  exphaz <- function(x,age,sex,year) { expectedhaz(ratetable, age=age, sex=sex, year=year, time=x)}
  
  survO <- as.matrix(predO$predictions)
  dimnames(survO) <- NULL
  
  N <- dim(survO)[1]
  P <- dim(survO)[2]
  
  hP <- matrix(-99, ncol=length(times), nrow=N)
  
  for (i in 1:N) # @Thomas : merci de voir si tu augmenter la vitesse du calcul de hP
  {
    hP[i,] <- sapply(times, FUN="exphaz", age=data[i,ratetable_vars$age],
                     sex=data[i,ratetable_vars$sex], year=data[i,ratetable_vars$year]) * inter
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
  
  # warning -> NA pour tCure ...
  
  res <- list(formula = formula,
              data = data,
              ratetable = ratetable,
              age = age,
              sex= sex,
              year = year,
              pro.time = pro.time,
              inter = splann$inter,
              size = splann$size,
              decay = splann$decay,
              fitsurvivalnet = splann,
              predictions = list(times=times, CIFc=distE, CIFp = distP, tPcure = Pcure, aPcure = distPinf, Sp=survP, Sc =  survU)
  )
  
  class(res) <- "rsPLANN"
  return(res)
}
