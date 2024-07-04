
print.survivalPLANN <- function (x, ...)
{    
   nmiss <- sum(x$missing)

   print(x$fitsurvivalnet, ...)
   
   cat("\n", length(x$intervals)-1 ," intervals computed, ranging from ", 
   min(x$intervals), " to ", max(x$interval), " with a step of ", x$inter, 
   " unit.", sep="")
   
   if(nmiss==1) { cat("\n", "(", nmiss, 
                      " observation deleted due to missingness)", sep="") }
   if(nmiss >1) { cat("\n", "(", nmiss,
                      " observations deleted due to missingness)", sep="") }
}



