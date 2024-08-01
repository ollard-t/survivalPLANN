
print.survivalPLANN <- function (x, ...)
{    
   nmiss <- sum(x$missing)
   
   message("\n", "A neural network based on the nnet package and the PLANN method ", sep="")
   
   message("\n", length(x$intervals)-1 ," intervals computed, ranging from ", 
   min(x$intervals), " to ", max(x$interval), " with a step of ", x$inter, 
   " unit.", sep="")
   
   if(nmiss==1) { message("\n", "(", nmiss, 
                      " observation deleted due to missingness)", sep="") }
   if(nmiss >1) { message("\n", "(", nmiss,
                      " observations deleted due to missingness)", sep="") }
}
